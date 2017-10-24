#include "PoissonRealBarcodesMergeStrategy.h"

#include <Tools/UtilFunctions.h>
#include <Tools/Logs.h>

namespace Estimation
{
namespace Merge
{
	PoissonRealBarcodesMergeStrategy::PoissonRealBarcodesMergeStrategy(const PoissonTargetEstimator &target_estimator,
	                                                                   const barcodes_parser_ptr barcodes_parser,
	                                                                   size_t min_genes_before_merge,
	                                                                   size_t min_genes_after_merge,
	                                                                   unsigned max_merge_edit_distance)
		: RealBarcodesMergeStrategy(barcodes_parser, min_genes_before_merge, min_genes_after_merge,
		                            max_merge_edit_distance, 0)
		, _target_estimator(target_estimator)
	{}

	unsigned PoissonRealBarcodesMergeStrategy::get_max_merge_dist(unsigned min_real_cb_dist) const
	{
		return min_real_cb_dist == 0 ? 2 : min_real_cb_dist + 1;
	}

	std::string PoissonRealBarcodesMergeStrategy::merge_type() const
	{
		return "Poisson";
	}

	void PoissonRealBarcodesMergeStrategy::init(const Estimation::CellsDataContainer &container)
	{
		RealBarcodesMergeStrategy::init(container);
		this->_target_estimator.init(container.umi_distribution());
	}

	void PoissonRealBarcodesMergeStrategy::release()
	{
		this->_target_estimator.release();
		RealBarcodesMergeStrategy::release();
	}

	long PoissonRealBarcodesMergeStrategy::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
																 const MergeStrategyAbstract::ul_list_t &neighbour_cells)
	{
		return this->_target_estimator.get_best_merge_target(container, base_cell_ind, neighbour_cells);
	}

	size_t PoissonRealBarcodesMergeStrategy::get_log_period() const
	{
		return 1000;
	}
}
}
