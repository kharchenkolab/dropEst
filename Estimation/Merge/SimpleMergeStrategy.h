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
		typedef boost::unordered_map<size_t, size_t> u_u_hash_t;
		typedef boost::unordered_set<size_t> sul_set_t;
		typedef boost::unordered_map<std::string, sul_set_t> sul_l_map_t;
		static const double EPS;

	private:
		size_t get_cells_with_common_umigs(Estimation::CellsDataContainer &container, const ul_list_t &cb_reassigned,
										   const Tools::IndexedValue &processed_genes_count,
										   const sul_l_map_t &umig_cell_ids, u_u_hash_t &umig_top) const;

		static sul_l_map_t get_umigs_map(const CellsDataContainer& container);

	protected:
		virtual void merge_inited(Estimation::CellsDataContainer &container, ul_list_t &filtered_cells) const override;

	public:
		SimpleMergeStrategy(int min_genes_before_merge, int min_genes_after_merge, int max_merge_edit_distance,
						  double min_merge_fraction);

		bool merge(Estimation::CellsDataContainer &container, long target_cell_ind, double target_cell_fraction,
				   const Tools::IndexedValue &source_genes_count, ul_list_t &cb_reassigned, ISIHM &cb_reassigned_to) const;

		virtual std::string merge_type() const override;
	};
}
}