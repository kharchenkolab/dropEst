#include "BamProcessor.h"

#include "Tools/Logs.h"

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>

namespace Estimation
{
	BamProcessor::BamProcessor(size_t read_prefix_length, const std::string &reads_params_name, const std::string &gtf_path)
		: _use_names_map(reads_params_name.length() != 0)
		, _read_prefix_length(read_prefix_length)
		, _use_genes_container(gtf_path.length() != 0)
	{
		if (this->_use_names_map)
		{
			std::ifstream ifs(reads_params_name);
			boost::archive::binary_iarchive ia(ifs);
			ia >> this->_reads_params;
		}

		if (this->_use_genes_container)
		{
			this->_genes_container = Tools::RefGenesContainer(gtf_path);
		}
	}

	void BamProcessor::parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
									   CellsDataContainer &container) const
	{
		CellsDataContainer::s_i_hash_t cb_ids;
		CellsDataContainer::s_ii_hash_t umig_cells_counts;

		for (auto const &file : bam_files)
		{
			this->parse_bam_file(file, print_result_bams, cb_ids, umig_cells_counts, container);
		}

		container.init_filled(umig_cells_counts);
	}

	void BamProcessor::parse_bam_file(const std::string &bam_name, bool print_result_bam, CellsDataContainer::s_i_hash_t &cells_ids,
							   CellsDataContainer::s_ii_hash_t &umig_cells_counts, CellsDataContainer &container) const
	{
		using namespace BamTools;
		L_TRACE << "Start reading bam files";

		BamReader reader;
		BamWriter writer;
		if (!reader.Open(bam_name))
			throw std::runtime_error("Could not open BAM file: " + bam_name);

		if (print_result_bam)
		{
			std::string result_bam_name = BamProcessor::get_result_bam_name(bam_name);
			if (!writer.Open(result_bam_name, reader.GetHeader(), reader.GetReferenceData()))
				throw std::runtime_error("Could not open BAM file to write: " + result_bam_name);
		}

		BamAlignment alignment;
		int max_cell_id = 0;
		long total_reads = 0, exonic_reads = 0;
		while (reader.GetNextAlignment(alignment))
		{
			total_reads++;
			if (total_reads % 1000000 == 0)
			{
				L_TRACE << "Total " << total_reads << " reads processed";
			}

			if (alignment.Length < this->_read_prefix_length)
			{
				L_ERR << "WARNING: read is shorter than read_prefix_length. total_reads=" << total_reads;
				continue;
			}

			std::string chr_name = reader.GetReferenceData()[alignment.RefID].RefName;

			Tools::ReadParameters read_params;
			const std::string read_name = alignment.Name;
			if (!this->get_read_params(read_name, read_params))
				continue;

			std::string cell_barcode = read_params.cell_barcode(), umi = read_params.umi_barcode();

			std::string gene = this->get_gene(chr_name, alignment);
			if (print_result_bam)
			{
				BamProcessor::write_alignment(writer, alignment, gene, read_params);
			}

			if (gene == "")
			{
				container.stats().inc_cell_chr_umi(chr_name, cell_barcode, Stats::NON_EXONE);
				continue;
			}

			exonic_reads++;
			L_DEBUG << read_name << " cell:" << cell_barcode << " UMI:" << umi << " prefix:"
					<< alignment.QueryBases.substr(this->_read_prefix_length) << "\tXF:" << gene;

			max_cell_id = std::max(BamProcessor::save_read_data(chr_name, cell_barcode, umi, gene, cells_ids,
																umig_cells_counts, container), max_cell_id);
		}

		reader.Close();
		writer.Close();
		L_TRACE << "Done (" << total_reads << " total reads; " << exonic_reads << " exonic reads; "
				<< max_cell_id + 1 << " cell barcodes)";
	}

	bool BamProcessor::get_read_params(const std::string &read_name, Tools::ReadParameters &read_params) const
	{
		if (this->_use_names_map)
		{
			auto iter = this->_reads_params.find(read_name);
			if (iter == this->_reads_params.end())
			{
				L_ERR << "WARNING: can't find read_name in map: " << read_name;
				return false;
			}

			read_params = iter->second;
			if (read_params.is_empty())
			{
				L_ERR << "WARNING: empty parameters for read_name: " << read_name;
				return false;
			}
		}
		else
		{
			try
			{
				read_params = Tools::ReadParameters(read_name);
			}
			catch (std::runtime_error &error)
			{
				L_ERR << error.what();
				return false;
			}
		}

		return true;
	}

	int BamProcessor::save_read_data(const std::string &chr_name, const std::string &cell_barcode, const std::string &umi,
									  const std::string &gene, CellsDataContainer::s_i_hash_t &cells_ids,
									  CellsDataContainer::s_ii_hash_t &umig_cells_counts, CellsDataContainer &container)
	{
		int cell_id = container.add_record(cell_barcode, umi, gene, cells_ids);

		std::string umig = umi + gene; // +iseq
		umig_cells_counts[umig][cell_id]++;

		L_DEBUG << "UMIg=" << umig_cells_counts[umig].size();

		container.stats().inc(Stats::READS_BY_UMIG, cell_barcode + "_" + umig);
		container.stats().inc(Stats::READS_BY_CB, cell_barcode);
		container.stats().inc_cell_chr_umi(chr_name, cell_barcode, Stats::EXONE);
		return cell_id;
	}

	std::string BamProcessor::get_result_bam_name(const std::string &bam_name)
	{
		return bam_name.substr(0, bam_name.find_last_of(".")) + ".corrected.bam";
	}

	std::string BamProcessor::get_gene(const std::string &chr_name, BamTools::BamAlignment alignment) const
	{
		if (this->_use_genes_container)
			return this->_genes_container.get_gene_info(chr_name, alignment.Position,
												 alignment.Position + alignment.Length).id();

		std::string gene;
		alignment.GetTag("GE", gene);
		return gene;
	}

	void BamProcessor::write_alignment(BamTools::BamWriter &writer, BamTools::BamAlignment &alignment,
									   const std::string &gene, const Tools::ReadParameters &parameters) const
	{
		alignment.Name = parameters.read_name();
		if (gene != "")
		{
			alignment.AddTag("GE", "Z", gene);
		}
		writer.SaveAlignment(alignment);
	}
}
