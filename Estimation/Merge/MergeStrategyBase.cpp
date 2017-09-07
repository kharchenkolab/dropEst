#include <Tools/Logs.h>
#include <Estimation/CellsDataContainer.h>
#include <Tools/IndexedValue.h>
#include "MergeStrategyBase.h"

#include <numeric>

namespace Estimation
{
namespace Merge
{
MergeStrategyAbstract::ul_list_t MergeStrategyBase::merge_inited(CellsDataContainer &container) const
{
	id_id_set_map_t cb_reassigned_to_it;
	ul_list_t cb_reassign_targets(container.total_cells_number());
	std::iota(cb_reassign_targets.begin(), cb_reassign_targets.end(), 0);

	auto const& filtered_cells = container.filtered_cells();
	std::vector<long> target_cell_inds(filtered_cells.size());

	for (size_t genes_count_id = 0; genes_count_id < filtered_cells.size(); ++genes_count_id)
	{
		target_cell_inds[genes_count_id] = this->get_merge_target(container, filtered_cells[genes_count_id]);
		if (genes_count_id % this->get_log_period() == 0 && genes_count_id > 0)
		{
			L_TRACE << "Total " << genes_count_id << " tags processed";
		}
	}

	size_t merges_count = 0, excluded_cells_num = 0;
	for (size_t genes_count_id = 0; genes_count_id < filtered_cells.size(); ++genes_count_id)
	{
		size_t base_cell_ind = filtered_cells[genes_count_id];
		long target_cell_ind = target_cell_inds[genes_count_id];
		if (target_cell_ind < 0)
		{
			container.exclude_cell(base_cell_ind);
			excluded_cells_num++;
			continue;
		}

		if (target_cell_ind != cb_reassign_targets.at(target_cell_ind))
		{
			target_cell_ind = cb_reassign_targets[target_cell_ind]; // For the case when real barcodes could be merged too
		}

		if (target_cell_ind == base_cell_ind)
			continue;

		this->merge_force(container, base_cell_ind, (size_t)target_cell_ind, cb_reassign_targets, cb_reassigned_to_it);
		merges_count++;
	}

	L_TRACE << "Total " << merges_count << " cells merged";
	L_TRACE << "Total " << excluded_cells_num << " cells excluded";

	return cb_reassign_targets;
}

	size_t MergeStrategyBase::get_log_period() const
	{
		return 100000;
	}

	void MergeStrategyBase::reassign(size_t cell_id, size_t target_cell_id, ul_list_t &cb_reassign_targets,
	                                 id_id_set_map_t &cb_reassigned_to_it) const
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

void MergeStrategyBase::merge_force(Estimation::CellsDataContainer &container, size_t src_cell_id,
                                    size_t target_cell_ind, ul_list_t &cb_reassign_targets,
                                    id_id_set_map_t &cb_reassigned_to_it) const
{
	container.merge_cells(src_cell_id, target_cell_ind);
	this->reassign(src_cell_id, target_cell_ind, cb_reassign_targets, cb_reassigned_to_it);
}

MergeStrategyBase::MergeStrategyBase(size_t min_genes_before_merge, size_t min_genes_after_merge,
                                     unsigned max_merge_edit_distance, double min_merge_fraction)
	: MergeStrategyAbstract(min_genes_before_merge, min_genes_after_merge)
	, _max_merge_edit_distance(max_merge_edit_distance)
	, _min_merge_fraction(min_merge_fraction)
{}


size_t MergeStrategyBase::get_umigs_intersect_size(const Cell &cell1, const Cell &cell2)
{
	std::map<std::string, Gene>::const_iterator gene1_it = cell1.genes().begin(); //Not unordered!!!
	std::map<std::string, Gene>::const_iterator gene2_it = cell2.genes().begin();

	size_t intersect_size = 0;
	while (gene1_it != cell1.genes().end() && gene2_it != cell2.genes().end())
	{
		int comp_res = gene1_it->first.compare(gene2_it->first);
		if (comp_res < 0)
		{
			gene1_it++;
			continue;
		}

		if (comp_res > 0)
		{
			gene2_it++;
			continue;
		}

		auto umi1_it = gene1_it->second.umis().begin();
		auto umi2_it = gene2_it->second.umis().begin();

		while (umi1_it != gene1_it->second.umis().end() && umi2_it != gene2_it->second.umis().end())
		{
			comp_res = umi1_it->first.compare(umi2_it->first);
			if (comp_res < 0)
			{
				++umi1_it;
				continue;
			}

			if (comp_res > 0)
			{
				++umi2_it;
				continue;
			}

			++intersect_size;
			++umi2_it;
			++umi1_it;
		}

		++gene1_it;
		++gene2_it;
	}

	return intersect_size;
}
}
}