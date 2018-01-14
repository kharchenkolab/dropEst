#pragma once

#include "MergeStrategyBase.h"

#include <Estimation/CellsDataContainer.h>
#include <Tools/UtilFunctions.h>

namespace Estimation
{
namespace Merge
{
	class SimpleMergeStrategy : public MergeStrategyBase
	{
	private:
		using sul_set_t = std::unordered_set<size_t>;
		using umig_map_t = std::unordered_map<std::pair<StringIndexer::index_t, StringIndexer::index_t>, sul_set_t, Tools::PairHash>;
		static const double EPS;

	private:
		umig_map_t _cell_ids_by_umig;

	protected:
		u_u_hash_t get_cells_with_common_umigs(const CellsDataContainer &container, size_t base_cell_ind) const;

		long get_merge_target(CellsDataContainer &container, size_t base_cell_ind) override;
		void init(const CellsDataContainer &container) override;
		void release() override;

	public:
		SimpleMergeStrategy(size_t min_genes_before_merge, size_t min_genes_after_merge,
		                    unsigned max_merge_edit_distance, double min_merge_fraction);

		std::string merge_type() const override;
	};
}
}