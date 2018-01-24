#pragma once

#include "MergeStrategyBase.h"

#include <Estimation/CellsDataContainer.h>
#include <Tools/UtilFunctions.h>
#include <limits>

namespace Estimation
{
namespace Merge
{
	class MergeAllMergeStrategy : public MergeStrategyBase
	{
	protected:
		long get_merge_target(CellsDataContainer &container, size_t base_cell_ind) override
		{
			int min_ed = std::numeric_limits<int>::max();
			int max_umi_num = 0;
			size_t target_ind = std::numeric_limits<size_t>::max();
			for (auto const &cell_ind: container.filtered_cells())
			{
				const size_t target_umi_num = container.cell(cell_ind).umis_number();
				if (target_umi_num <= container.cell(base_cell_ind).umis_number())
					continue;

				int ed = Tools::edit_distance(container.cell(base_cell_ind).barcode_c(),
				                              container.cell(cell_ind).barcode_c(), false, this->_max_merge_edit_distance);

				if (ed > this->_max_merge_edit_distance)
					continue;

				if (min_ed > ed)
				{
					min_ed = ed;
					max_umi_num = target_umi_num;
					target_ind = cell_ind;
				}
				else if (min_ed == ed & max_umi_num < target_umi_num)
				{
					max_umi_num = target_umi_num;
					target_ind = cell_ind;
				}
			}

			if (target_ind != std::numeric_limits<size_t>::max())
				return target_ind;

			return base_cell_ind;
		}

	public:
		MergeAllMergeStrategy(size_t min_genes_before_merge, size_t min_genes_after_merge,
		                      unsigned max_merge_edit_distance)
			: MergeStrategyBase(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, 0)
		{}

		std::string merge_type() const override
		{
			return "Merge all";
		}
	};
}
}
