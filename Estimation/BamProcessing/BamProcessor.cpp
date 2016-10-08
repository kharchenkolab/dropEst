#include <Tools/ReadParameters.h>
#include <Tools/Logs.h>
#include <Estimation/Stats.h>
#include <iomanip>
#include "BamProcessor.h"
#include "BamController.h"

namespace Estimation
{
	namespace BamProcessing
	{
		BamProcessor::BamProcessor(CellsDataContainer &container, bool print_bam)
				: container(container)
				, print_bam(print_bam)
				, total_exonic_reads(0)
				, total_reads(0)
		{}

		BamProcessor::~BamProcessor()
		{
			this->writer.Close();
		}

		void BamProcessor::save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi, const std::string& gene)
		{
			if (gene == "")
			{
				this->container.stats().inc(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, cell_barcode, chr_name);
				return;
			}

			size_t cell_id = this->container.add_record(cell_barcode, umi, gene, this->cells_ids);
			if (this->container.cell_genes(cell_id).at(gene).at(umi) == 1)
			{
				this->container.stats().inc(Stats::EXONE_UMIS_PER_CHR_PER_CELL, cell_barcode, chr_name);
			}

			this->container.stats().inc(Stats::EXONE_READS_PER_CHR_PER_CELL, cell_barcode, chr_name);
			this->total_exonic_reads++;
		}

		void BamProcessor::inc_reads()
		{
			this->total_reads++;
		}

		void BamProcessor::trace_state(const std::string &bam_file) const
		{
			L_TRACE << bam_file << ": " << this->total_reads << " total reads; " << this->total_exonic_reads << std::setprecision(3)
					<< " ("<< (100.0*this->total_exonic_reads / this->total_reads) <<"%) exonic; " << this->cells_ids.size() << " CBs";
		}

		void BamProcessor::update_bam(const std::string &bam_file, const BamTools::BamReader &reader)
		{
			if (!this->print_bam)
				return;

			this->writer.Close();
			std::string result_bam_name = this->get_result_bam_name(bam_file);
			if (!this->writer.Open(result_bam_name, reader.GetHeader(), reader.GetReferenceData()))
				throw std::runtime_error("Could not open BAM file to write: " + result_bam_name);
		}

		std::string BamProcessor::get_result_bam_name(const std::string &bam_name) const
		{
			return bam_name.substr(0, bam_name.find_last_of(".")) + ".tagged.bam";
		}

		void BamProcessor::write_alignment(BamTools::BamAlignment alignment, const std::string &gene,
										const Tools::ReadParameters &read_params)
		{
			if (!this->print_bam)
				return;

			alignment.Name = read_params.read_name_safe();
			if (gene != "")
			{
				alignment.AddTag(BamController::GENE_TAG, "Z", gene);
			}

			alignment.AddTag(BamController::CB_TAG, "Z", read_params.cell_barcode());
			alignment.AddTag(BamController::UMI_TAG, "Z", read_params.umi_barcode());
			this->writer.SaveAlignment(alignment);
		}
	}
}