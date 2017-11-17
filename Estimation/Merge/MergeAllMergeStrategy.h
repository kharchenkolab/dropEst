#pragma once

#include "MergeStrategyBase.h"

#include <Estimation/CellsDataContainer.h>
#include <Tools/UtilFunctions.h>

namespace Estimation
{
namespace Merge
{
	class MergeAllMergeStrategy : public MergeStrategyBase
	{
	protected:
		virtual long get_merge_target(CellsDataContainer &container, size_t base_cell_ind) override
		{
			for (auto const &cell_ind: container.filtered_cells())
			{
				if (container.cell(cell_ind).umis_number() <= container.cell(base_cell_ind).umis_number())
					continue;

				int ed = Tools::edit_distance(container.cell(base_cell_ind).barcode_c(),
				                              container.cell(cell_ind).barcode_c(), false, this->_max_merge_edit_distance);

				if (ed > this->_max_merge_edit_distance)
					continue;

				return cell_ind;
			}

			return base_cell_ind;
		}

	public:
		MergeAllMergeStrategy(size_t min_genes_before_merge, size_t min_genes_after_merge,
		                      unsigned max_merge_edit_distance)
			: MergeStrategyBase(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, 0)
		{}

		virtual std::string merge_type() const override
		{
			return "Merge all";
		}
	};
}
}