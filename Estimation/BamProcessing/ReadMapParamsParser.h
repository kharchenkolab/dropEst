#pragma once
#include "BamProcessor.h"
#include "ReadsParamsParser.h"
#include <Tools/ReadParameters.h>
#include <api/BamAlignment.h>

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadMapParamsParser : public ReadsParamsParser
		{
		private:
			Tools::reads_params_map_t _reads_params;

		private:
			void init(const std::string & read_param_names, bool save_read_name);

		public:
			ReadMapParamsParser(const std::string &gtf_path, bool save_read_names,
								const std::string &read_param_names, int gene_match_level);
			virtual ~ReadMapParamsParser();

			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) override;
		};
	}
}