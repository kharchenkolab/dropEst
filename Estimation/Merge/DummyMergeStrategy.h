#pragma once

#include "MergeStrategyAbstract.h"

namespace Estimation
{
namespace Merge
{
	class DummyMergeStrategy : public MergeStrategyAbstract
	{
	protected:
		virtual ul_list_t merge_inited(Estimation::CellsDataContainer &container) override
		{
			ul_list_t reassign_targets(container.total_cells_number());
			std::iota(reassign_targets.begin(), reassign_targets.end(), 0);
			return reassign_targets;
		}

	public:
		DummyMergeStrategy(size_t min_genes_before_merge, size_t min_genes_after_merge)
			: MergeStrategyAbstract(min_genes_before_merge, min_genes_after_merge)
		{}

		virtual std::string merge_type() const override {
			return "No";
		}
	};
}
}