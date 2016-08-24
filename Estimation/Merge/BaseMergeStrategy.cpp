#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>
#include "BaseMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	const double BaseMergeStrategy::EPS = 0.00001;

	using Tools::IndexedValue;

	void BaseMergeStrategy::merge(Estimation::CellsDataContainer &container, const s_ii_hash_t &umig_cells_counts,
								  ids_t &filtered_cells) const
	{
		int merges_count = 0;

		ISIHM cb_reassigned_to_it;
		ids_t cb_reassign_targets(container.cell_barcodes().size());
		std::iota(cb_reassign_targets.begin(), cb_reassign_targets.end(), 0);

		L_TRACE << "merging linked tags ";

		int tag_index = 0;
		for (auto const &genes_count : container.cells_genes_counts_sorted())
		{ // iterate through the minimally-selected CBs, from low to high counts
			if (++tag_index % 1000 == 0)
			{
				L_TRACE << "Total " << tag_index << " tags processed, " << merges_count << " cells merged";
			}

			i_i_hash_t umigs_intersect_top;
			size_t umigs_count = this->get_umigs_intersect_top(container, cb_reassign_targets, genes_count,
															   umig_cells_counts, umigs_intersect_top);

			int top_cell_ind = -1;
			double top_cb_fraction = -1;
			long top_cb_genes_count = -1;
			for (auto const &cell: umigs_intersect_top)
			{
				int cell_ind = cell.first;
				double cb_fraction = cell.second / (double) umigs_count;
				if (cb_fraction - top_cb_fraction > EPS || (abs(cb_fraction - top_cb_fraction) < EPS &&
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
				container.stats().add_merge_count(genes_count.value);
				if (genes_count.value >= this->min_genes_after_merge())
				{
					filtered_cells.push_back(genes_count.index);
				}
			}
		}

		container.stats().merge(cb_reassign_targets, container.cell_barcodes());

		if (filtered_cells.size() > 1)
		{
			std::reverse(filtered_cells.begin(), filtered_cells.end());
		}

		L_INFO << "Done (" << merges_count << " merges performed)" << std::endl;
	}

	bool BaseMergeStrategy::merge(Estimation::CellsDataContainer &container, size_t target_cell_ind,
								  double target_cell_fraction, const IndexedValue &source_genes_count,
								  ids_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it) const
	{
		if (target_cell_ind < 0)
			return false;

		// check if the top candidate is valid for merging
		if (target_cell_fraction < this->_min_merge_fraction)
			return false;

		size_t source_cell_ind = source_genes_count.index;
		int ed = Tools::edit_distance(container.cell_barcode(target_cell_ind).c_str(),
									  container.cell_barcode(source_cell_ind).c_str());
		if (ed >= this->_max_merge_edit_distance)
			return false;

		this->merge_force(container, source_cell_ind, (size_t)target_cell_ind, source_genes_count.value,
						  cb_reassign_targets, cb_reassigned_to_it);

		return true;
	}

	size_t BaseMergeStrategy::get_umigs_intersect_top(Estimation::CellsDataContainer &container,
													  const ids_t &cb_reassign_targets, const IndexedValue &processed_genes_count,
													  const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const
	{
		size_t umigs_count = 0;
		for (auto const &gene: container.cell_genes(processed_genes_count.index))
		{
			const std::string &gene_name = gene.first;
			const s_i_map_t &umis = gene.second;

			for (auto const &umi_count: umis)
			{
				std::string umig = umi_count.first + gene_name;
				const i_i_hash_t &umig_cells = umigs_cells_counts.at(umig);
				for (auto const &cell : umig_cells)
				{
					int cell_with_same_umig_id = cell.first;
					size_t reassign_target = cb_reassign_targets[cell_with_same_umig_id]; //if not reassigned then cell_with_same_umig_id
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

	BaseMergeStrategy::BaseMergeStrategy(int min_genes_before_merge, int min_genes_after_merge,
										 int max_merge_edit_distance, double min_merge_fraction)
			: AbstractMergeStrategy(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance,
									min_merge_fraction)
	{}
}
}