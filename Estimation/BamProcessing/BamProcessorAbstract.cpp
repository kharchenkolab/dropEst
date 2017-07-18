#include <Tools/Logs.h>
#include "BamProcessorAbstract.h"
#include "BamController.h"

namespace Estimation
{
namespace BamProcessing
{
	BamProcessorAbstract::BamProcessorAbstract()
		: _total_reads_num(0)
		, _cant_parse_reads_num(0)
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

	void BamProcessorAbstract::save_alignment(BamTools::BamAlignment alignment, const std::string &name,
											  const std::string &gene, const std::string &barcode, const std::string &umi)
	{
		alignment.Name = name;
		if (gene != "")
		{
			alignment.AddTag(BamController::GENE_TAG, "Z", gene);
		}

		alignment.AddTag(BamController::CB_TAG, "Z", barcode);
		alignment.AddTag(BamController::UMI_TAG, "Z", umi);
		this->_writer.SaveAlignment(alignment);
	}

	void BamProcessorAbstract::inc_cant_parse_num()
	{
		++this->_cant_parse_reads_num;
	}
}
}
