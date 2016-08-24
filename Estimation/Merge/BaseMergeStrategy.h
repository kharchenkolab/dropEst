#pragma once

#include "AbstractMergeStrategy.h"

#include <Estimation/CellsDataContainer.h>

namespace Estimation
{
namespace Merge
{
	class BaseMergeStrategy : public AbstractMergeStrategy
	{
	private:
		static const double EPS;

	private:
		size_t get_umigs_intersect_top(Estimation::CellsDataContainer &container, const ids_t &cb_reassigned,
									   const Tools::IndexedValue &processed_genes_count,
									   const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const;

	public:
		BaseMergeStrategy(int min_genes_before_merge, int min_genes_after_merge, int max_merge_edit_distance,
						  double min_merge_fraction);

		virtual void merge(Estimation::CellsDataContainer &container, const s_ii_hash_t &umig_cells_counts,
						   ids_t &filtered_cells) const override;

		bool merge(Estimation::CellsDataContainer &container, size_t target_cell_ind, double target_cell_fraction,
				   const Tools::IndexedValue &source_genes_count, ids_t &cb_reassigned, ISIHM &cb_reassigned_to) const;
	};
}
}