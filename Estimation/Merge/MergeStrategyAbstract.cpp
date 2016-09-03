#include <Tools/Logs.h>
#include "MergeStrategyAbstract.h"

namespace Estimation
{
namespace Merge
{

	MergeStrategyAbstract::MergeStrategyAbstract(int min_genes_before_merge, int min_genes_after_merge)
			: _min_genes_before_merge(min_genes_before_merge)
			, _min_genes_after_merge(min_genes_after_merge)
	{}

	void MergeStrategyAbstract::merge(Estimation::CellsDataContainer &container, const s_uu_hash_t &umig_cells_counts,
									  ul_list_t &filtered_cells)
	{
		Tools::trace_time("Start merge");
		this->init(container);
		Tools::trace_time("Merge initialized");
		this->merge_inited(container, umig_cells_counts, filtered_cells);
		this->release();
		Tools::trace_time("Merge finished");
	}

	void MergeStrategyAbstract::init(const Estimation::CellsDataContainer &container)
	{}

	void MergeStrategyAbstract::release()
	{}

	int MergeStrategyAbstract::min_genes_before_merge() const
	{
		return this->_min_genes_before_merge;
	}

	int MergeStrategyAbstract::min_genes_after_merge() const
	{
		return this->_min_genes_after_merge;
	}
}
}