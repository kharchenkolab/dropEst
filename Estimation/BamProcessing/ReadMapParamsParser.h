#pragma once
#include "BamProcessor.h"
#include "ReadParamsParser.h"
#include <Tools/ReadParameters.h>
#include <api/BamAlignment.h>

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadMapParamsParser : public ReadParamsParser
		{
		private:
			Tools::reads_params_map_t _reads_params;

		private:
			void init(const std::string &read_param_filenames);

		public:
			ReadMapParamsParser(const std::string &gtf_path, const std::string &read_param_filenames,
			                    const BamTags &tags, bool gene_in_chromosome_name);
			virtual ~ReadMapParamsParser();

			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) override;
		};
	}
}