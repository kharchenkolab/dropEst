#include "BamProcessorFactory.h"

#include "Estimation/BamProcessing/BamProcessor.h"
#include "Estimation/BamProcessing/FilledBamProcessor.h"
#include "Estimation/BamProcessing/ReadMapBamProcessor.h"

namespace Estimation
{
namespace BamProcessing
{
	std::shared_ptr<BamProcessor> BamProcessorFactory::get(bool filled_bam, const std::string &reads_params_names_str,
														   const std::string &gtf_filename, size_t read_prefix_length)
	{
		if (filled_bam)
			return std::make_shared<FilledBamProcessor>(read_prefix_length, gtf_filename);

		if (reads_params_names_str != "")
			return std::make_shared<ReadMapBamProcessor>(read_prefix_length, reads_params_names_str, gtf_filename);

		return std::make_shared<BamProcessor>(read_prefix_length, gtf_filename);
	}
}
}