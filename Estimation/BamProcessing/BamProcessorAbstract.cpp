#include "BamProcessorAbstract.h"

#include "BamController.h"
#include <Tools/Logs.h>
#include <Tools/ReadParameters.h>

namespace Estimation
{
namespace BamProcessing
{
	BamProcessorAbstract::BamProcessorAbstract(const BamTags &tags_info)
		: _total_reads_num(0)
		, _cant_parse_reads_num(0)
		, _low_quality_reads_num(0)
		, _tags(tags_info)
	{}

	BamProcessorAbstract::~BamProcessorAbstract()
	{
		this->_writer.Close();
	}

	void BamProcessorAbstract::inc_reads()
	{
		this->_total_reads_num++;
	}

	void BamProcessorAbstract::update_bam(const std::string &bam_file, const BamTools::BamReader &reader)
	{
		this->_writer.Close();
		std::string result_bam_name = this->get_result_bam_name(bam_file);
		std::string::size_type path_end = result_bam_name.find_last_of("\\/");
		if (path_end != std::string::npos)
		{
			result_bam_name = result_bam_name.substr(path_end + 1);
		}

		L_TRACE << "Write to " << result_bam_name;

		if (!this->_writer.Open(result_bam_name, reader.GetHeader(), reader.GetReferenceData()))
			throw std::runtime_error("Could not open BAM file to write: " + result_bam_name);
	}

	size_t BamProcessorAbstract::total_reads_num() const
	{
		return this->_total_reads_num;
	}

	size_t BamProcessorAbstract::cant_parse_reads_num() const
	{
		return this->_cant_parse_reads_num;
	}

	size_t BamProcessorAbstract::low_quality_reads_num() const
	{
		return this->_low_quality_reads_num;
	}

	void BamProcessorAbstract::save_alignment(BamTools::BamAlignment alignment, const ReadInfo &read_info_raw,
	                                          const std::string &cell_barcode_corrected, const std::string &umi_corrected)
	{
		auto const &raw_params = read_info_raw.params;

		if (!read_info_raw.gene.empty())
		{
			alignment.EditTag(this->_tags.gene, "Z", read_info_raw.gene);
		}

		alignment.EditTag(this->_tags.cell_barcode_raw, "Z", raw_params.cell_barcode());
		alignment.EditTag(this->_tags.umi_raw, "Z", raw_params.umi());

		if (!raw_params.cell_barcode_quality().empty())
		{
			alignment.EditTag(this->_tags.cell_barcode_quality, "Z", raw_params.cell_barcode_quality());
		}

		if (!raw_params.umi_quality().empty())
		{
			alignment.EditTag(this->_tags.umi_quality, "Z", raw_params.umi_quality());
		}

		// Read type
		if (read_info_raw.umi_mark == UMI::Mark::HAS_EXONS)
		{
			alignment.EditTag(this->_tags.read_type, "Z", this->_tags.exonic_read_value_out);
		}
		else if (read_info_raw.umi_mark == UMI::Mark::HAS_INTRONS)
		{
			alignment.EditTag(this->_tags.read_type, "Z", this->_tags.intronic_read_value_out);
		}
		else if (read_info_raw.umi_mark == UMI::Mark::HAS_NOT_ANNOTATED)
		{
			alignment.EditTag(this->_tags.read_type, "Z", this->_tags.intergenic_read_value_out);
		}

		// Corrected barcodes
		if (!cell_barcode_corrected.empty())
		{
			alignment.EditTag(this->_tags.cell_barcode, "Z", cell_barcode_corrected);
		}

		if (!umi_corrected.empty())
		{
			alignment.EditTag(this->_tags.umi, "Z", umi_corrected);
		}

		this->_writer.SaveAlignment(alignment);
	}

	void BamProcessorAbstract::inc_cant_parse_num()
	{
		++this->_cant_parse_reads_num;
	}

	void BamProcessorAbstract::inc_low_quality_num()
	{
		++this->_low_quality_reads_num;
	}
}
}
