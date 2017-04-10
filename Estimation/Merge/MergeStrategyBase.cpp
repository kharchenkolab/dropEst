#include <Tools/Logs.h>
#include <Estimation/CellsDataContainer.h>
#include <Tools/IndexedValue.h>
#include <boost/range/adaptor/reversed.hpp>
#include "MergeStrategyBase.h"

namespace Estimation
{
namespace Merge
{
MergeStrategyAbstract::ul_list_t MergeStrategyBase::merge_inited(CellsDataContainer &container) const
{
	ISIHM cb_reassigned_to_it;
	ul_list_t cb_reassign_targets(container.cell_barcodes_raw().size());
	std::iota(cb_reassign_targets.begin(), cb_reassign_targets.end(), 0);

	auto const& cell_genes_counts = container.cells_gene_counts_sorted();
	std::vector<long> target_cell_inds(cell_genes_counts.size());

	for (size_t genes_count_id = 0; genes_count_id < cell_genes_counts.size(); ++genes_count_id)
	{
		target_cell_inds[genes_count_id] = this->get_merge_target(container, cell_genes_counts[genes_count_id].index);
		if (genes_count_id % this->get_log_period() == 0 && genes_count_id > 0)
		{
			L_TRACE << "Total " << genes_count_id << " tags processed";
		}
	}

	size_t merges_count = 0;
	for (size_t genes_count_id = 0; genes_count_id < cell_genes_counts.size(); ++genes_count_id)
	{
		size_t base_cell_ind = cell_genes_counts[genes_count_id].index;
		long target_cell_ind = target_cell_inds[genes_count_id];
		if (target_cell_ind < 0)
		{
			container.exclude_cell(base_cell_ind);
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

	size_t excluded_cells_num = container.excluded_cells().size();
	L_INFO << "Total " << merges_count << " merges";
	L_INFO << "Total " << excluded_cells_num << " cells excluded";
	L_INFO << container.cells_gene_counts_sorted().size() - merges_count - excluded_cells_num << " cells with " << this->min_genes_before_merge() << " genes left";

	container.update_cell_sizes(this->min_genes_after_merge(), false);

	return cb_reassign_targets;
}

	size_t MergeStrategyBase::get_log_period() const
	{
		return 100000;
	}

	void MergeStrategyBase::reassign(size_t cell_id, size_t target_cell_id, ul_list_t &cb_reassign_targets,
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

void MergeStrategyBase::merge_force(Estimation::CellsDataContainer &container, size_t src_cell_id,
									    size_t target_cell_ind, ul_list_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it) const
{
	L_DEBUG << "Merge: " << container.cell_barcode(src_cell_id) << " to " << container.cell_barcode(target_cell_ind);
	container.stats().inc(Estimation::Stats::MERGES_COUNT_PER_CB, container.cell_barcode(target_cell_ind));
	container.stats().set(Estimation::Stats::MERGE_TARGET_BY_BASE, container.cell_barcode(src_cell_id), container.cell_barcode(target_cell_ind), 1);

	container.merge_cells(src_cell_id, target_cell_ind);
	this->reassign(src_cell_id, target_cell_ind, cb_reassign_targets, cb_reassigned_to_it);
}

MergeStrategyBase::MergeStrategyBase(const boost::property_tree::ptree &config)
	: MergeStrategyAbstract(config)
	, _max_merge_edit_distance(config.get<unsigned>("max_merge_edit_distance"))
	, _min_merge_fraction(config.get<double>("min_merge_fraction"))
{}


size_t MergeStrategyBase::get_umigs_intersect_size(const genes_t &cell1_dist, const genes_t &cell2_dist)
{
	std::map<std::string, CellsDataContainer::s_i_map_t>::const_iterator gene1_it = cell1_dist.begin(); //Not unordered!!!
	std::map<std::string, CellsDataContainer::s_i_map_t>::const_iterator gene2_it = cell2_dist.begin();

	size_t intersect_size = 0;
	while (gene1_it != cell1_dist.end() && gene2_it != cell2_dist.end())
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

		auto umi1_it = gene1_it->second.begin();
		auto umi2_it = gene2_it->second.begin();

		while (umi1_it != gene1_it->second.end() && umi2_it != gene2_it->second.end())
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