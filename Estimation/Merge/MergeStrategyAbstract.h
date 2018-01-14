#pragma once

#include <Estimation/CellsDataContainer.h>
#include "Tools/IndexedValue.h"

namespace Estimation
{
class CellsDataContainer;

namespace Merge
{
	class MergeStrategyAbstract
	{
	public:
		using ul_list_t = std::vector<size_t>;

	private:
		const size_t _min_genes_before_merge;
		const size_t _min_genes_after_merge;

	protected:
		virtual ul_list_t merge_inited(Estimation::CellsDataContainer &container) = 0;

	public:
		MergeStrategyAbstract(size_t min_genes_before_merge, size_t min_genes_after_merge);

		ul_list_t merge(Estimation::CellsDataContainer &container);

		virtual void init(const Estimation::CellsDataContainer &container);
		virtual void release();
		virtual std::string merge_type() const = 0;

		size_t min_genes_before_merge() const;
		size_t min_genes_after_merge() const;
	};
}
}
