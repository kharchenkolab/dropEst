#include <Tools/UtilFunctions.h>
#include "PoissonSimpleMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
PoissonSimpleMergeStrategy::PoissonSimpleMergeStrategy(int min_genes_before_merge,
													   int min_genes_after_merge,
													   int max_merge_edit_distance,
													   double min_merge_fraction)
		: SimpleMergeStrategy(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
{}

long PoissonSimpleMergeStrategy::get_merge_target(const Estimation::CellsDataContainer &container, size_t base_cell_ind) const
{
	u_u_hash_t cells_with_common_umigs = this->get_cells_with_common_umigs(container, base_cell_ind);

	ul_list_t neighbour_cells;
	neighbour_cells.reserve(cells_with_common_umigs.size());
	for (auto const &cell : cells_with_common_umigs)
	{
		if (Tools::edit_distance(container.cell_barcode(base_cell_ind).c_str(),
								 container.cell_barcode(cell.first).c_str()) > this->_max_merge_edit_distance)
			continue;

		neighbour_cells.push_back(cell.first);
	}

	if (neighbour_cells.size() == 0)
		return base_cell_ind;

	return this->_target_estimator.get_best_merge_target(container, base_cell_ind, neighbour_cells);
}

void PoissonSimpleMergeStrategy::init(const CellsDataContainer &container)
{
	SimpleMergeStrategy::init(container);
	this->_target_estimator.init(container);
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
}
}