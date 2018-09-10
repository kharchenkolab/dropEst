#pragma once
#include "BamProcessor.h"
#include "ReadParamsParser.h"
#include <Estimation/StringIndexer.h>
#include <Estimation/ReadParametersEfficient.h>
#include <Tools/ReadParameters.h>
#include <api/BamAlignment.h>

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadMapParamsParser : public ReadParamsParser
		{
		private:
			char _min_barcode_phred;
			ReadParametersEfficient::reads_params_map_t _read_params;

			StringIndexer _barcode_indexer;
			StringIndexer _umi_indexer;
			StringIndexer _umi_quality_indexer;

		private:
			void init(const std::string &read_param_filenames);

		public:
			ReadMapParamsParser(const std::string &gtf_path, const std::string &read_param_filenames,
			                    const BamTags &tags, bool gene_in_chromosome_name, char min_barcode_phred);

			bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) override;
		};
	}
}