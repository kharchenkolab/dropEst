#include "MergeStrategyFactory.h"

#include "FilteringMergeStrategy.h"
#include "SimpleMergeStrategy.h"
#include "RealBarcodesMergeStrategy.h"
#include "PoissonRealBarcodesMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	std::shared_ptr<MergeStrategyAbstract> MergeStrategyFactory::get(double min_merge_fraction, unsigned min_genes_before_merge,
											   unsigned min_genes_after_merge, unsigned max_merge_edit_distance,
											   bool merge_tags, const std::string &barcodes_filename,
											   size_t barcode2_length, bool use_old_merge)
	{
		if (!merge_tags)
			return std::shared_ptr<MergeStrategyAbstract>(new FilteringMergeStrategy(min_genes_before_merge, min_genes_after_merge));

		if (barcodes_filename == "")
			return std::shared_ptr<MergeStrategyAbstract>(new SimpleMergeStrategy(min_genes_before_merge, min_genes_after_merge,
																		 max_merge_edit_distance, min_merge_fraction));
		if (use_old_merge)
			return std::shared_ptr<MergeStrategyAbstract>(new RealBarcodesMergeStrategy(barcodes_filename, barcode2_length,
																			 min_genes_before_merge, min_genes_after_merge,
																			 max_merge_edit_distance, min_merge_fraction));

		return std::shared_ptr<MergeStrategyAbstract>(new PoissonRealBarcodesMergeStrategy(barcodes_filename, barcode2_length,
																					min_genes_before_merge, min_genes_after_merge,
																					max_merge_edit_distance, min_merge_fraction));
	}
}
}