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
		typedef boost::unordered_set<size_t> sul_set_t;
		typedef boost::unordered_map<std::pair<StringIndexer::index_t, StringIndexer::index_t>, sul_set_t> umig_map_t;
		static const double EPS;

	private:
		umig_map_t _umig_cell_ids;

	protected:
		u_u_hash_t get_cells_with_common_umigs(const CellsDataContainer &container, size_t base_cell_ind) const;

		virtual long get_merge_target(const CellsDataContainer &container, size_t base_cell_ind) const override;
		virtual void init(const CellsDataContainer &container) override;
		virtual void release() override;

	public:
		SimpleMergeStrategy(size_t min_genes_before_merge, size_t min_genes_after_merge,
		                    unsigned max_merge_edit_distance, double min_merge_fraction);

		virtual std::string merge_type() const override;
	};
}
}