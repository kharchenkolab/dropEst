#include "MergeStrategyFactory.h"

#include "FilteringMergeStrategy.h"
#include "BaseMergeStrategy.h"
#include "RealBarcodesMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	std::shared_ptr<IMergeStrategy> MergeStrategyFactory::get(double min_merge_fraction, unsigned min_genes_before_merge,
											   unsigned min_genes_after_merge, unsigned max_merge_edit_distance,
											   bool merge_tags, const std::string &barcodes_filename,
											   size_t barcode2_length)
	{
		if (!merge_tags)
			return std::shared_ptr<IMergeStrategy>(new FilteringMergeStrategy(min_genes_before_merge, min_genes_after_merge));

		if (barcodes_filename == "")
			return std::shared_ptr<IMergeStrategy>(new BaseMergeStrategy(min_genes_before_merge, min_genes_after_merge,
																		 max_merge_edit_distance, min_merge_fraction));

		return std::shared_ptr<IMergeStrategy>(new RealBarcodesMergeStrategy(barcodes_filename, barcode2_length,
																			 min_genes_before_merge, min_genes_after_merge,
																			 max_merge_edit_distance, min_merge_fraction));
	}
}
}