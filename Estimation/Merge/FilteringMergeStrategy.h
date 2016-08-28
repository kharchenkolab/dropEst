#pragma once

#include "MergeStrategyAbstract.h"

#include <boost/range/adaptor/reversed.hpp>

namespace Estimation
{
namespace Merge
{
	class FilteringMergeStrategy : public MergeStrategyAbstract
	{
	protected:
		virtual void merge_inited(Estimation::CellsDataContainer &container, const s_uu_hash_t &umig_cells_counts,
							  ul_list_t &filtered_cells) const override
		{
			for (auto const &gene_count : boost::adaptors::reverse(container.cells_genes_counts_sorted()))
			{
				if (gene_count.value < this->min_genes_after_merge())
					break;

				filtered_cells.push_back(gene_count.index);
			}
		}

	public:
		FilteringMergeStrategy(unsigned min_genes_before_merge, unsigned min_genes_after_merge)
				: MergeStrategyAbstract(min_genes_before_merge, min_genes_after_merge)
		{}
	};
}
}