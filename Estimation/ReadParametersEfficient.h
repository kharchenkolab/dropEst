#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <Tools/ReadParameters.h>
#include "StringIndexer.h"

namespace Estimation
{
	class ReadParametersEfficient
	{
	public:
		using reads_params_map_t = std::unordered_map<std::string, ReadParametersEfficient>;

	private:
		size_t _cell_barcode_id;
		size_t _umi_id;
		size_t _umi_quality_id;
		bool _pass_quality_threshold;

		bool _is_empty;

	public:
		explicit ReadParametersEfficient(const Tools::ReadParameters &parameters, StringIndexer &barcode_indexer,
		                                 StringIndexer &umi_indexer, StringIndexer &umi_quality_indexer);

		Tools::ReadParameters parameters(StringIndexer &barcode_indexer, StringIndexer &umi_indexer,
		                                 StringIndexer &umi_quality_indexer) const;
	};
}
