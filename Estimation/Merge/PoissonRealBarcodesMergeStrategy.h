#pragma once

#include <RInside.h>
#include "RealBarcodesMergeStrategy.h"
#include "PoissonTargetEstimator.h"

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeRejections;
}

namespace Estimation
{
namespace Merge
{
	class PoissonRealBarcodesMergeStrategy : public RealBarcodesMergeStrategy
	{
		friend struct TestEstimatorMergeProbs::testPoissonMergeRejections;

	private:
		PoissonTargetEstimator _target_estimator;

	protected:
		void init(const Estimation::CellsDataContainer &container) override;
		void release() override;

		long get_best_merge_target(CellsDataContainer &container, size_t base_cell_ind,
		                                   const ul_list_t &neighbour_cells) override;
		unsigned get_max_merge_dist(unsigned min_real_cb_dist) const override;

		size_t get_log_period() const override;

	public:
		PoissonRealBarcodesMergeStrategy(const PoissonTargetEstimator &target_estimator,
		                                 const barcodes_parser_ptr barcodes_parser, size_t min_genes_before_merge,
		                                 size_t min_genes_after_merge, unsigned max_merge_edit_distance);

		std::string merge_type() const override;
	};
}
}