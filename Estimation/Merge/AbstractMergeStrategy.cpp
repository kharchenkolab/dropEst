#include <Tools/Logs.h>
#include <Estimation/CellsDataContainer.h>
#include "AbstractMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
void AbstractMergeStrategy::reassign(size_t cell_id, size_t target_cell_id, ids_t &cb_reassign_targets,
											ISIHM &cb_reassigned_to_it) const
{
	cb_reassign_targets[cell_id] = target_cell_id; // set reassignment mapping
	cb_reassigned_to_it[target_cell_id].insert(cell_id); // reassign current cell

	// transfer mapping of the cbs previously mapped to kid
	auto reassigned_to_cell_iter = cb_reassigned_to_it.find(cell_id);
	if (reassigned_to_cell_iter == cb_reassigned_to_it.end())
		return;

	for (auto reassigned_id: reassigned_to_cell_iter->second)
	{
		cb_reassign_targets[reassigned_id] = target_cell_id; // update reassignment mapping
		cb_reassigned_to_it[target_cell_id].insert(reassigned_id);
	}

	reassigned_to_cell_iter->second.clear();
}

void AbstractMergeStrategy::merge_force(Estimation::CellsDataContainer &container, size_t src_cell_id,
									    size_t target_cell_ind, ids_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it) const
{
	L_DEBUG << "Merge: " << container.cell_barcode(src_cell_id) << " to " << container.cell_barcode(target_cell_ind);
	container.stats().inc(Estimation::Stats::MERGES_COUNT_PER_CB, container.cell_barcode(target_cell_ind));

	container.merge(src_cell_id, target_cell_ind);
	this->reassign(src_cell_id, target_cell_ind, cb_reassign_targets, cb_reassigned_to_it);
}

AbstractMergeStrategy::AbstractMergeStrategy(int min_genes_before_merge, int min_genes_after_merge,
													int max_merge_edit_distance, double min_merge_fraction)
		: IMergeStrategy(min_genes_before_merge, min_genes_after_merge)
		, _max_merge_edit_distance(max_merge_edit_distance)
		, _min_merge_fraction(min_merge_fraction)
{}

}
}