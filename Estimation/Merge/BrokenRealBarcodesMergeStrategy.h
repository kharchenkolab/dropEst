#pragma once

#include <RInside.h>
#include "RealBarcodesMergeStrategy.h"

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeProbs;
	struct testPoissonMergeInit;
	struct testPoissonMergeTime;
}

namespace Estimation
{
namespace Merge
{
	class [[deprecated("Was used for test purposes only")]] BrokenRealBarcodesMergeStrategy : public RealBarcodesMergeStrategy
	{
	protected:
		virtual void init(const Estimation::CellsDataContainer &container) override;

		virtual long get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind, const ul_list_t &neighbour_cells) const override;
		virtual unsigned get_max_merge_dist(unsigned min_real_cb_dist) const override;

	public:
		BrokenRealBarcodesMergeStrategy(barcodes_parser_ptr barcodes_parser,
		                                unsigned min_genes_before_merge, unsigned min_genes_after_merge,
		                                unsigned max_merge_edit_distance, double min_merge_fraction);

		virtual std::string merge_type() const override;
	};
}
}