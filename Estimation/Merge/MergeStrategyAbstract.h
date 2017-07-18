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
		typedef std::vector<size_t> ul_list_t;

	private:
		const int _min_genes_before_merge;
		const int _min_genes_after_merge;

	protected:
		virtual ul_list_t merge_inited(Estimation::CellsDataContainer &container) const = 0;

	public:
		MergeStrategyAbstract(unsigned min_genes_before_merge, unsigned min_genes_after_merge);

		ul_list_t merge(Estimation::CellsDataContainer &container);

		virtual void init(const Estimation::CellsDataContainer &container);
		virtual void release();
		virtual std::string merge_type() const = 0;

		int min_genes_before_merge() const;
		int min_genes_after_merge() const;
	};
}
}
