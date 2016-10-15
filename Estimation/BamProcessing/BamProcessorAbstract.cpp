#include "BamProcessorAbstract.h"

namespace Estimation
{
namespace BamProcessing
{
	BamProcessorAbstract::BamProcessorAbstract()
			: _total_reads(0)
	{}

	BamProcessorAbstract::~BamProcessorAbstract()
	{
		this->_writer.Close();
	}

	void BamProcessorAbstract::inc_reads()
	{
		this->_total_reads++;
	}

	void BamProcessorAbstract::update_bam(const std::string &bam_file, const BamTools::BamReader &reader)
	{
		this->_writer.Close();
		std::string result_bam_name = this->get_result_bam_name(bam_file);
		if (!this->_writer.Open(result_bam_name, reader.GetHeader(), reader.GetReferenceData()))
			throw std::runtime_error("Could not open BAM file to write: " + result_bam_name);
	}

	void BamProcessorAbstract::save_alignment(BamTools::BamAlignment alignment)
	{
		this->_writer.SaveAlignment(alignment);
	}

	size_t BamProcessorAbstract::total_reads() const
	{
		return this->_total_reads;
	}
}
}
