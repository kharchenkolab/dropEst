#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>
#include "SimpleMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	const double SimpleMergeStrategy::EPS = 0.00001;

	using Tools::IndexedValue;

	void SimpleMergeStrategy::merge_inited(Estimation::CellsDataContainer &container, ul_list_t &filtered_cells) const
	{
		sul_l_map_t umig_cell_ids = SimpleMergeStrategy::get_umigs_map(container);

		size_t merges_count = 0;

		ISIHM cb_reassigned_to_it;
		ul_list_t cb_reassign_targets(container.cell_barcodes_raw().size());
		std::iota(cb_reassign_targets.begin(), cb_reassign_targets.end(), 0);

		L_TRACE << "merging linked tags ";

		int tag_index = 0;
		for (auto const &genes_count : container.cells_genes_counts_sorted())
		{ // iterate through the minimally-selected CBs, from low to high counts
			if (++tag_index % 1000 == 0)
			{
				L_TRACE << "Total " << tag_index << " tags processed, " << merges_count << " cells merged";
			}

			u_u_hash_t cells_with_common_umigs;
			size_t umigs_count = this->get_cells_with_common_umigs(container, cb_reassign_targets, genes_count,
																   umig_cell_ids, cells_with_common_umigs);

			long top_cell_ind = -1;
			double top_cb_fraction = -1;
			long top_cb_genes_count = -1;
			for (auto const &cell: cells_with_common_umigs)
			{
				size_t cell_ind = cell.first;
				double cb_fraction = cell.second / (double) umigs_count;
				if (cb_fraction - top_cb_fraction > EPS || (std::abs(cb_fraction - top_cb_fraction) < EPS &&
															container.cell_genes(cell_ind).size() > top_cb_genes_count))
				{
					top_cell_ind = cell_ind;
					top_cb_fraction = cb_fraction;
					top_cb_genes_count = container.cell_genes(cell_ind).size();
				}
			}

			if (this->merge(container, top_cell_ind, top_cb_fraction, genes_count, cb_reassign_targets,
							cb_reassigned_to_it))
			{
				merges_count++;
			}
			else
			{
				if (genes_count.value >= this->min_genes_after_merge())
				{
					filtered_cells.push_back(genes_count.index);
				}
			}
		}

		container.stats().merge(cb_reassign_targets, container.cell_barcodes_raw());

		if (filtered_cells.size() > 1)
		{
			std::reverse(filtered_cells.begin(), filtered_cells.end());
		}

		L_INFO << "Done (" << merges_count << " merges performed)" << std::endl;
	}

	bool SimpleMergeStrategy::merge(Estimation::CellsDataContainer &container, long target_cell_ind,
								  double target_cell_fraction, const IndexedValue &source_genes_count,
								  ul_list_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it) const
	{
		if (target_cell_ind < 0)
			return false;

		// check if the top candidate is valid for merging
		if (target_cell_fraction < this->_min_merge_fraction)
			return false;

		size_t source_cell_ind = source_genes_count.index;
		int ed = Tools::edit_distance(container.cell_barcode((size_t)target_cell_ind).c_str(),
									  container.cell_barcode(source_cell_ind).c_str());
		if (ed >= this->_max_merge_edit_distance)
			return false;

		this->merge_force(container, source_cell_ind, (size_t)target_cell_ind, cb_reassign_targets, cb_reassigned_to_it);

		return true;
	}

	size_t SimpleMergeStrategy::get_cells_with_common_umigs(Estimation::CellsDataContainer &container,
															const ul_list_t &cb_reassign_targets,
															const IndexedValue &processed_genes_count,
															const sul_l_map_t &umig_cell_ids, u_u_hash_t &umig_top) const
	{
		size_t umigs_count = 0;
		for (auto const &gene: container.cell_genes(processed_genes_count.index))
		{
			const std::string &gene_name = gene.first;
			const s_i_map_t &umis = gene.second;

			for (auto const &umi_count: umis)
			{
				std::string umig = umi_count.first + gene_name;
				const auto &umig_cells = umig_cell_ids.at(umig);
				for (size_t cell_with_same_umig_id : umig_cells)
				{
					size_t reassign_target = cb_reassign_targets.at(cell_with_same_umig_id); //if not reassigned then cell_with_same_umig_id
					if (reassign_target == processed_genes_count.index)
						continue;

					if (container.cell_genes(reassign_target).size() > processed_genes_count.value)
					{
						umig_top[reassign_target]++;
					}
				}
				umigs_count++;
			}
		}
		return umigs_count;
	}

	SimpleMergeStrategy::SimpleMergeStrategy(int min_genes_before_merge, int min_genes_after_merge,
										 int max_merge_edit_distance, double min_merge_fraction)
			: MergeStrategyBase(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance,
									min_merge_fraction)
	{}

	std::string SimpleMergeStrategy::merge_type() const
	{
		return "Simple";
	}

	SimpleMergeStrategy::sul_l_map_t SimpleMergeStrategy::get_umigs_map(const CellsDataContainer& container)
	{
		sul_l_map_t umig_cell_ids;
		for (auto const &genes_count : container.cells_genes_counts_sorted())
		{
			for (auto const &gene : container.cell_genes(genes_count.index))
			{
				for (auto const &umi : gene.second)
				{
					auto res = umig_cell_ids.emplace(std::make_pair(umi.first + gene.first, sul_set_t()));
					res.first->second.emplace(genes_count.index);
				}
			}
		}

		return umig_cell_ids;
	}
}
}