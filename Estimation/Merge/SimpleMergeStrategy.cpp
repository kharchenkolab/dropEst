#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>
#include "SimpleMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	const double SimpleMergeStrategy::EPS = 0.00001;

	SimpleMergeStrategy::SimpleMergeStrategy(size_t min_genes_before_merge, size_t min_genes_after_merge,
	                                         unsigned max_merge_edit_distance, double min_merge_fraction)
			: MergeStrategyBase(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
	{}

	SimpleMergeStrategy::u_u_hash_t SimpleMergeStrategy::get_cells_with_common_umigs(const CellsDataContainer &container,
																					 size_t base_cell_ind) const
	{
		u_u_hash_t common_umigs_per_cell;
		for (auto const &gene: container.cell(base_cell_ind).genes())
		{
			for (auto const &umi_count: gene.second.umis())
			{
				auto umig = std::make_pair(umi_count.first, gene.first);
				const auto &umig_cells = this->_umig_cell_ids.at(umig);
				for (size_t cell_with_same_umig_id : umig_cells)
				{
					if (cell_with_same_umig_id == base_cell_ind)
						continue;

					// We don't need to check if cell is reassigned, because of sorting in descending order
					if (container.cell(cell_with_same_umig_id).size() >= container.cell(base_cell_ind).size())
					{
						common_umigs_per_cell[cell_with_same_umig_id]++;
					}
				}
			}
		}

		return common_umigs_per_cell;
	}

	std::string SimpleMergeStrategy::merge_type() const
	{
		return "Simple";
	}

	long SimpleMergeStrategy::get_merge_target(const CellsDataContainer &container, size_t base_cell_ind)
	{
		u_u_hash_t cells_with_common_umigs = this->get_cells_with_common_umigs(container, base_cell_ind);

		long top_cell_ind = -1;
		double top_cb_fraction = -1;
		long top_cb_genes_count = -1;
		for (auto const &cell: cells_with_common_umigs)
		{
			size_t cell_ind = cell.first;
			double cb_fraction =  0.5 * cell.second *
					(1. / container.cell(base_cell_ind).umis_number() + 1. / container.cell(cell_ind).umis_number());

//			container.stats().set(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
//								  container.cell_barcode(cell_ind), cb_fraction);

			if (cb_fraction - top_cb_fraction > SimpleMergeStrategy::EPS ||
					(std::abs(cb_fraction - top_cb_fraction) < SimpleMergeStrategy::EPS &&
							container.cell(cell_ind).size() > top_cb_genes_count))
			{
				int ed = Tools::edit_distance(container.cell(base_cell_ind).barcode_c(),
											  container.cell(cell_ind).barcode_c());

				if (ed >= this->_max_merge_edit_distance)
					continue;

//				container.stats().set(Stats::MERGE_EDIT_DISTANCE_BY_CELL, container.cell_barcode(base_cell_ind),
//									  container.cell_barcode(cell_ind), ed);

				top_cell_ind = cell_ind;
				top_cb_fraction = cb_fraction;
				top_cb_genes_count = container.cell(cell_ind).size();
			}
		}

		if (top_cb_fraction < this->_min_merge_fraction)
			return base_cell_ind;

		return top_cell_ind;
	}

	void SimpleMergeStrategy::init(const CellsDataContainer &container)
	{
		MergeStrategyAbstract::init(container);
		for (auto const &cell_id : container.filtered_cells())
		{
			for (auto const &gene : container.cell(cell_id).genes())
			{
				for (auto const &umi : gene.second.umis())
				{
					auto umig = std::make_pair(umi.first, gene.first);
					auto res = this->_umig_cell_ids.emplace(std::make_pair(umig, sul_set_t()));
					res.first->second.emplace(cell_id);
				}
			}
		}
	}

	void SimpleMergeStrategy::release()
	{
		this->_umig_cell_ids.clear();
		MergeStrategyAbstract::release();
	}
}
}
