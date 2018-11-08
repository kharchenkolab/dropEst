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
		{}

		void BamProcessor::save_read(const ReadInfo &read_info)
		{
			this->_container.add_record(read_info);
		}

		void BamProcessor::trace_state(const std::string &trace_prefix) const
		{
			std::stringstream cant_parse_msg;
			if (this->cant_parse_reads_num() > 0)
			{
				cant_parse_msg << ", can't parse " << (100.0*this->cant_parse_reads_num() / this->total_reads_num()) <<"% reads, ";
			}

			if (this->low_quality_reads_num() > 0)
			{
				cant_parse_msg << "low-quality reads: " << (100.0*this->low_quality_reads_num() / this->total_reads_num()) <<"%, ";
			}

			if (this->nonmapped_reads_num() > 0)
			{
				cant_parse_msg << "low-quality reads: " << (100.0*this->low_quality_reads_num() / this->total_reads_num()) <<"%";
			}

			L_TRACE << trace_prefix << ": " << this->total_reads_num() << " primarily aligned reads (" << std::setprecision(3)
					<< (100.0*this->container().intergenic_reads_num() / this->total_reads_num()) <<"% intergenic, "
					<< (100.0*this->container().has_exon_reads_num() / this->total_reads_num()) <<"% touch exon, "
					<< (100.0*this->container().has_intron_reads_num() / this->total_reads_num()) <<"% touch intron, "
					<< (100.0*this->container().has_not_annotated_reads_num() / this->total_reads_num()) <<"% touch both gene and not annotated regions"
					<< cant_parse_msg.str()
					<< "), "
					<< this->nonmapped_reads_num() << " non-mapped reads, "
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

		void BamProcessor::write_alignment(BamTools::BamAlignment alignment, const ReadInfo &read_info)
		{
			if (!this->_print_bam)
				return;

			this->save_alignment(alignment, read_info);
		}

		const CellsDataContainer &BamProcessor::container() const
		{
			return this->_container;
		}
	}
}
