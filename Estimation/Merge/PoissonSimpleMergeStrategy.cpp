#include <Tools/UtilFunctions.h>
#include "PoissonSimpleMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
PoissonSimpleMergeStrategy::PoissonSimpleMergeStrategy(const PoissonTargetEstimator &target_estimator,
                                                       unsigned min_genes_before_merge, unsigned min_genes_after_merge,
                                                       unsigned max_merge_edit_distance)
	: SimpleMergeStrategy(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, 0)
	, _target_estimator(target_estimator)
{}

long PoissonSimpleMergeStrategy::get_merge_target(CellsDataContainer &container, size_t base_cell_ind)
{
	const std::string &base_cb = container.cell(base_cell_ind).barcode();
	u_u_hash_t cells_with_common_umigs = this->get_cells_with_common_umigs(container, base_cell_ind);

	ul_list_t neighbour_cells;
	neighbour_cells.reserve(cells_with_common_umigs.size());
	for (auto const &cell : cells_with_common_umigs)
	{
		const std::string &current_cb = container.cell(cell.first).barcode();
		unsigned edit_distance = Tools::edit_distance(base_cb.c_str(), current_cb.c_str());

		if (edit_distance > this->_max_merge_edit_distance)
			continue;

//		container.stats().set(Stats::MERGE_EDIT_DISTANCE_BY_CELL, base_cb, current_cb, edit_distance);
		neighbour_cells.push_back(cell.first);
	}

	if (neighbour_cells.size() == 0)
		return base_cell_ind;

	long target = this->_target_estimator.get_best_merge_target(container, base_cell_ind, neighbour_cells);
	if (target != -1)
		return target;

	return base_cell_ind;
}

void PoissonSimpleMergeStrategy::init(const CellsDataContainer &container)
{
	SimpleMergeStrategy::init(container);
	this->_target_estimator.init(container.umi_distribution());
}

void PoissonSimpleMergeStrategy::release()
{
	this->_target_estimator.release();
	SimpleMergeStrategy::release();
}

std::string PoissonSimpleMergeStrategy::merge_type() const
{
	return "Poisson Simple";
}

size_t PoissonSimpleMergeStrategy::get_log_period() const
{
	return 10000;
}
}
}
