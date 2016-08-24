#pragma once

#include "IMergeStrategy.h"

#include <boost/range/adaptor/reversed.hpp>

namespace Estimation
{
namespace Merge
{
	class FilteringMergeStrategy : public IMergeStrategy
	{
	public:
		FilteringMergeStrategy(unsigned min_genes_before_merge, unsigned min_genes_after_merge)
				: IMergeStrategy(min_genes_before_merge, min_genes_after_merge)
		{}

		virtual void merge(Estimation::CellsDataContainer &container, const s_ii_hash_t &umig_cells_counts,
						   ids_t &filtered_cells) const override
		{
			for (auto const &gene_count : boost::adaptors::reverse(container.cells_genes_counts_sorted()))
			{
				if (gene_count.value < this->min_genes_after_merge())
					break;

				filtered_cells.push_back(gene_count.index);
			}
		}
	};
}
}