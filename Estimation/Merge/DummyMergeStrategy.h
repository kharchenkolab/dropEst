#pragma once

#include "MergeStrategyAbstract.h"

#include <boost/range/adaptor/reversed.hpp>

namespace Estimation
{
namespace Merge
{
	class DummyMergeStrategy : public MergeStrategyAbstract
	{
	protected:
		virtual ul_list_t merge_inited(Estimation::CellsDataContainer &container) const override
		{
			ul_list_t reassign_targets(container.cell_barcodes_raw().size());
			std::iota(reassign_targets.begin(), reassign_targets.end(), 0);
			return reassign_targets;
		}

	public:
		DummyMergeStrategy(unsigned min_genes_before_merge, unsigned min_genes_after_merge)
				: MergeStrategyAbstract(min_genes_before_merge, min_genes_after_merge)
		{}

		virtual std::string merge_type() const override {
			return "No";
		}
	};
}
}