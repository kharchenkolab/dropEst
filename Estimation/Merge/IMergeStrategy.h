#pragma once

#include <Estimation/CellsDataContainer.h>
#include "Tools/IndexedValue.h"

namespace Estimation
{
class CellsDataContainer;

namespace Merge
{
	class IMergeStrategy
	{
	public:
		typedef boost::unordered_map<int, int> i_i_hash_t;
		typedef boost::unordered_map<std::string, i_i_hash_t> s_ii_hash_t;

		typedef std::map<std::string, int> s_i_map_t;

		typedef std::vector<Tools::IndexedValue> i_counter_t;
		typedef std::vector<long> ints_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<std::string> names_t;

	private:
		const int _min_genes_before_merge;
		const int _min_genes_after_merge;

	public:
		IMergeStrategy(int min_genes_before_merge, int min_genes_after_merge)
				: _min_genes_before_merge(min_genes_before_merge)
				, _min_genes_after_merge(min_genes_after_merge)
		{}

		virtual void merge(Estimation::CellsDataContainer &container, const s_ii_hash_t &umig_cells_counts,
						   ids_t &filtered_cells) const = 0;

		int min_genes_before_merge() const
		{
			return this->_min_genes_before_merge;
		}

		int min_genes_after_merge() const
		{
			return this->_min_genes_after_merge;
		}
	};
}
}