#pragma once

#include "IMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	class MergeStrategyFactory
	{
	public:
		static std::shared_ptr<IMergeStrategy> get(double min_merge_fraction, unsigned min_genes_before_merge,
												   unsigned min_genes_after_merge, unsigned max_merge_edit_distance,
												   bool merge_tags, const std::string &barcodes_filename,
												   size_t barcode2_length);
	};
}
}