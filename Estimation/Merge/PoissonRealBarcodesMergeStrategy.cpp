#include "PoissonRealBarcodesMergeStrategy.h"

#include <Tools/UtilFunctions.h>
#include <Tools/Logs.h>

namespace Estimation
{
namespace Merge
{
	PoissonRealBarcodesMergeStrategy::PoissonRealBarcodesMergeStrategy(const std::string &barcodes_filename,
																	   const boost::property_tree::ptree &config)
		: RealBarcodesMergeStrategy(barcodes_filename, config)
	{ }

	long PoissonRealBarcodesMergeStrategy::get_max_merge_dist(long min_real_cb_dist) const
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
		this->_target_estimator.init(container);
	}

	void PoissonRealBarcodesMergeStrategy::release()
	{
		this->_target_estimator.release();
		RealBarcodesMergeStrategy::release();
	}

	long PoissonRealBarcodesMergeStrategy::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
																 const MergeStrategyAbstract::ul_list_t &neighbour_cells) const
	{
		return this->_target_estimator.get_best_merge_target(container, base_cell_ind, neighbour_cells);
	}
}
}
