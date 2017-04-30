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
				: _container(container)
				, print_bam(print_bam)
				, total_intragenic_reads(0)
		{}

		void BamProcessor::save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi,
		                             const std::string& gene, const CellsDataContainer::Mark &umi_mark)
		{
			if (gene == "")
			{
				this->_container.stats().inc(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, cell_barcode, chr_name);
				this->total_intragenic_reads++;
				return;
			}

			size_t cell_id = this->_container.add_record(cell_barcode, umi, gene, umi_mark);
			if (this->_container.cell_genes(cell_id).at(gene).at(umi).read_count == 1)
			{
				this->_container.stats().inc(Stats::EXONE_UMIS_PER_CHR_PER_CELL, cell_barcode, chr_name);
			}

			this->_container.stats().inc(Stats::EXONE_READS_PER_CHR_PER_CELL, cell_barcode, chr_name);
		}

		void BamProcessor::trace_state(const std::string &trace_prefix) const
		{
			L_TRACE << trace_prefix << ": " << this->total_reads() << " total reads; " << std::setprecision(3)
					<< (100.0*this->total_intragenic_reads / this->total_reads()) <<"% intragenic; "
					<< (100.0*this->container().has_exon_reads_num() / this->total_reads()) <<"% touch exon; "
					<< (100.0*this->container().has_intron_reads_num() / this->total_reads()) <<"% touch intron; "
					<< (100.0*this->container().has_not_annotated_reads_num() / this->total_reads()) <<"% touch not annotated regions; "
					<< this->_container.cell_barcodes_raw().size() << " CBs read";
		}

		void BamProcessor::update_bam(const std::string &bam_file, const BamTools::BamReader &reader)
		{
			if (!this->print_bam)
				return;

			BamProcessorAbstract::update_bam(bam_file, reader);
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

			this->save_alignment(alignment, read_params.read_name_safe(), gene,
								 read_params.cell_barcode(), read_params.umi());
		}

		const CellsDataContainer &BamProcessor::container() const
		{
			return this->_container;
		}
	}
}