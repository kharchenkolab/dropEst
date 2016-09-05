#include "BamProcessor.h"

#include "Estimation/Stats.h"
#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"

#include <api/BamReader.h>
#include <api/BamWriter.h>

namespace Estimation
{
namespace BamProcessing
{
	const std::string BamProcessor::GENE_TAG = "GE";
	const std::string BamProcessor::CB_TAG = "CB";
	const std::string BamProcessor::UMI_TAG = "UB";

	BamProcessor::BamProcessor(size_t read_prefix_length, const std::string &gtf_path)
			: _read_prefix_length(read_prefix_length)
	{
		if (gtf_path.length() != 0)
		{
			this->_genes_container = Tools::RefGenesContainer(gtf_path);
		}
	}

	CellsDataContainer::s_uu_hash_t BamProcessor::parse_bam_files(const std::vector<std::string> &bam_files,
																  bool print_result_bams, CellsDataContainer &container) const
	{
		CellsDataContainer::s_i_map_t cb_ids;
		CellsDataContainer::s_uu_hash_t umig_cells_counts;

		for (auto const &file : bam_files)
		{
			this->parse_bam_file(file, print_result_bams, cb_ids, umig_cells_counts, container);
		}

		return umig_cells_counts;
	}

	void BamProcessor::parse_bam_file(const std::string &bam_name, bool print_result_bam, CellsDataContainer::s_i_map_t &cells_ids,
									  CellsDataContainer::s_uu_hash_t &umig_cells_counts, CellsDataContainer &container) const
	{
		this->init_temporaries_before_parsing(print_result_bam);

		using namespace BamTools;
		L_TRACE << "Start reading bam file: " + bam_name;

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
		std::unordered_set<std::string> unexpected_chromosomes;

		while (reader.GetNextAlignment(alignment))
		{
			total_reads++;
			if (total_reads % 1000000 == 0)
			{
				L_TRACE << "Total " << total_reads << " reads processed (" << exonic_reads << " exonic reads)";
			}

			if (alignment.Length < this->_read_prefix_length)
			{
				L_WARN << "WARNING: read is shorter than read_prefix_length. total_reads=" << total_reads;
				continue;
			}

			std::string chr_name = reader.GetReferenceData()[alignment.RefID].RefName;

			Tools::ReadParameters read_params;
			const std::string &read_name = alignment.Name;
			if (!this->get_read_params(alignment, read_params))
				continue;

			std::string cell_barcode = read_params.cell_barcode(), umi = read_params.umi_barcode();
			container.stats().inc(Stats::READS_PER_UMI_PER_CELL, cell_barcode, umi);

			std::string gene;
			try
			{
				gene = this->get_gene(chr_name, alignment);
			}
			catch (Tools::RefGenesContainer::ChrNotFoundException ex)
			{
				if (unexpected_chromosomes.emplace(ex.chr_name).second)
				{
					L_WARN << "WARNING: Can't find chromosome '" << ex.chr_name << "'";
				}
				continue;
			}

			if (print_result_bam)
			{
				BamProcessor::write_alignment(writer, alignment, gene, read_params);
			}

			if (gene == "")
			{
				L_DEBUG << "NonEx: " << read_name << ", cell: " << cell_barcode << " UMI: " << umi << ", start: " <<
						alignment.Position;
				container.stats().inc(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, cell_barcode, chr_name);
				continue;
			}

			exonic_reads++;
			L_DEBUG << "Exonic: " << read_name << ", cell: " << cell_barcode << " UMI: " << umi << ", prefix: "
					<< alignment.QueryBases.substr(this->_read_prefix_length) << "\tGene:" << gene;

			max_cell_id = std::max(BamProcessor::save_read_data(chr_name, cell_barcode, umi, gene, cells_ids,
																umig_cells_counts, container), max_cell_id);
		}

		reader.Close();
		writer.Close();

		this->release_temporaries_after_parsing();
		L_TRACE << "Done (" << total_reads << " total reads; " << exonic_reads << " exonic reads; "
				<< max_cell_id + 1 << " cell barcodes)";
	}

	bool BamProcessor::get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) const
	{
		try
		{
			read_params = Tools::ReadParameters(alignment.Name);
		}
		catch (std::runtime_error &error)
		{
			L_ERR << error.what();
			return false;
		}

		return true;
	}

	int BamProcessor::save_read_data(const std::string &chr_name, const std::string &cell_barcode, const std::string &umi,
									 const std::string &gene, CellsDataContainer::s_i_map_t &cells_ids,
									 CellsDataContainer::s_uu_hash_t &umig_cells_counts, CellsDataContainer &container)
	{
		int cell_id = container.add_record(cell_barcode, umi, gene, cells_ids);

		std::string umig = umi + gene;
		umig_cells_counts[umig][cell_id]++;

		L_DEBUG << "UMIg=" << umig_cells_counts[umig].size();

		container.stats().inc(Stats::EXONE_READS_PER_CB, cell_barcode);

		container.stats().inc(Stats::READS_PER_UMIG_PER_CELL, cell_barcode, umig);
		container.stats().inc(Stats::EXONE_READS_PER_CHR_PER_CELL, cell_barcode, chr_name);
		container.stats().inc(Stats::UMI_PER_CELL, cell_barcode, umi);
		return cell_id;
	}

	std::string BamProcessor::get_result_bam_name(const std::string &bam_name)
	{
		return bam_name.substr(0, bam_name.find_last_of(".")) + ".tagged.bam";
	}

	std::string BamProcessor::get_gene(const std::string &chr_name, BamTools::BamAlignment alignment) const
	{
		if (!this->_genes_container.is_empty())
			return this->_genes_container.get_gene_info(chr_name, alignment.Position,
														alignment.Position + alignment.Length).id();

		std::string gene;
		if (!alignment.GetTag(BamProcessor::GENE_TAG, gene))
			return "";

		return gene;
	}

	void BamProcessor::write_alignment(BamTools::BamWriter &writer, BamTools::BamAlignment &alignment,
									   const std::string &gene, const Tools::ReadParameters &parameters) const
	{
		alignment.Name = parameters.read_name_safe();
		if (gene != "")
		{
			alignment.AddTag(BamProcessor::GENE_TAG, "Z", gene);
		}

		alignment.AddTag(BamProcessor::CB_TAG, "Z", parameters.cell_barcode());
		alignment.AddTag(BamProcessor::UMI_TAG, "Z", parameters.umi_barcode());
		writer.SaveAlignment(alignment);
	}

	void BamProcessor::init_temporaries_before_parsing(bool save_read_name) const
	{}

	void BamProcessor::release_temporaries_after_parsing() const
	{}
}
}
