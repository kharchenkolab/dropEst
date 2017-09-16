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
		BamProcessor::BamProcessor(CellsDataContainer &container, const BamTags &tags, bool print_bam)
			: BamProcessorAbstract(tags)
			, _container(container)
			, _print_bam(print_bam)
			, _total_intergenic_reads(0)
		{}

		void BamProcessor::save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi,
		                             const std::string& gene, const UMI::Mark &umi_mark)
		{
			if (gene == "")
			{
				this->_total_intergenic_reads++;
			}

			this->_container.add_record(cell_barcode, umi, gene, chr_name, umi_mark);
		}

		void BamProcessor::trace_state(const std::string &trace_prefix) const
		{
			std::stringstream cant_parse_msg;
			if (this->cant_parse_reads_num() > 0)
			{
				cant_parse_msg << "can't parse " << (100.0*this->cant_parse_reads_num() / this->total_reads_num()) <<"% reads; ";
			}

			L_TRACE << trace_prefix << ": " << this->total_reads_num() << " total reads; " << std::setprecision(3)
					<< (100.0*this->_total_intergenic_reads / this->total_reads_num()) <<"% intergenic; "
					<< (100.0*this->container().has_exon_reads_num() / this->total_reads_num()) <<"% touch exon; "
					<< (100.0*this->container().has_intron_reads_num() / this->total_reads_num()) <<"% touch intron; "
					<< (100.0*this->container().has_not_annotated_reads_num() / this->total_reads_num()) <<"% touch not annotated regions; "
					<< cant_parse_msg.str()
					<< this->_container.total_cells_number() << " CBs read";
		}

		void BamProcessor::update_bam(const std::string &bam_file, const BamTools::BamReader &reader)
		{
			if (!this->_print_bam)
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
			if (!this->_print_bam)
				return;

			this->save_alignment(alignment, read_params, gene);
		}

		const CellsDataContainer &BamProcessor::container() const
		{
			return this->_container;
		}
	}
}
