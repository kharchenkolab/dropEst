#include <Tools/Logs.h>
#include "MergeStrategyAbstract.h"

namespace Estimation
{
namespace Merge
{
	MergeStrategyAbstract::MergeStrategyAbstract(unsigned min_genes_before_merge, unsigned min_genes_after_merge)
			: _min_genes_before_merge(min_genes_before_merge)
			, _min_genes_after_merge(std::max(min_genes_after_merge, min_genes_before_merge))
	{}

	MergeStrategyAbstract::ul_list_t MergeStrategyAbstract::merge(CellsDataContainer &container)
	{
		Tools::trace_time("Start merge");
		this->init(container);
		Tools::trace_time("Merge initialized");
		ul_list_t targets = this->merge_inited(container);
		this->release();
		Tools::trace_time("Merge finished");

		return targets;
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
