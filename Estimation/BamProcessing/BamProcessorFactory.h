#pragma once

#include <memory>
#include <string>

namespace Estimation
{
	namespace BamProcessing
	{
		class BamProcessor;
		class BamProcessorFactory
		{
		public:
			static std::shared_ptr<BamProcessor> get(bool filled_bam, const std::string &reads_params_names_str,
											  const std::string &gtf_filename, size_t read_prefix_length);
		};
	}
}