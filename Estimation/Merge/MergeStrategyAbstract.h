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
		typedef std::map<std::string, size_t> s_i_map_t;

		typedef std::vector<Tools::IndexedValue> i_counter_t;
		typedef std::vector<size_t> ul_list_t;
		typedef std::vector<std::string> names_t;

	private:
		const int _min_genes_before_merge;
		const int _min_genes_after_merge;

	protected:
		virtual void merge_inited(Estimation::CellsDataContainer &container, ul_list_t &filtered_cells) const = 0;

	public:
		MergeStrategyAbstract(int min_genes_before_merge, int min_genes_after_merge);

		void merge(Estimation::CellsDataContainer &container, ul_list_t &filtered_cells);

		virtual void init(const Estimation::CellsDataContainer &container);
		virtual void release();
		virtual std::string merge_type() const = 0;

		int min_genes_before_merge() const;
		int min_genes_after_merge() const;
	};
}
}