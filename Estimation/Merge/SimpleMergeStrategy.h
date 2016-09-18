#pragma once

#include "MergeStrategyBase.h"

#include <Estimation/CellsDataContainer.h>

namespace Estimation
{
namespace Merge
{
	class SimpleMergeStrategy : public MergeStrategyBase
	{
	private:
		static const double EPS;

	private:
		size_t get_umigs_intersect_top(Estimation::CellsDataContainer &container, const ul_list_t &cb_reassigned,
									   const Tools::IndexedValue &processed_genes_count,
									   const s_uu_hash_t &umigs_cells_counts, u_u_hash_t &umig_top) const;

	protected:
		virtual void merge_inited(Estimation::CellsDataContainer &container, const s_uu_hash_t &umig_cells_counts,
								  ul_list_t &filtered_cells) const override;

	public:
		SimpleMergeStrategy(int min_genes_before_merge, int min_genes_after_merge, int max_merge_edit_distance,
						  double min_merge_fraction);

		bool merge(Estimation::CellsDataContainer &container, long target_cell_ind, double target_cell_fraction,
				   const Tools::IndexedValue &source_genes_count, ul_list_t &cb_reassigned, ISIHM &cb_reassigned_to) const;

		virtual std::string merge_type() const override;
	};
}
}