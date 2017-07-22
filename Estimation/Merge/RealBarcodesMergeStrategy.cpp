#include "RealBarcodesMergeStrategy.h"

#include "BarcodesParsing/BarcodesParser.h"
#include "BarcodesParsing/ConstLengthBarcodesParser.h"
#include "BarcodesParsing/InDropBarcodesParser.h"

#include <fstream>

namespace Estimation
{
	namespace Merge
	{
		RealBarcodesMergeStrategy::RealBarcodesMergeStrategy(barcodes_parser_ptr barcodes_parser,
		                                                     unsigned min_genes_before_merge, unsigned min_genes_after_merge,
		                                                     unsigned max_merge_edit_distance, double min_merge_fraction)
				: MergeStrategyBase(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
				, _barcodes_parser(barcodes_parser)
		{
			this->_barcodes_parser->init();
		}

		long RealBarcodesMergeStrategy::get_merge_target(const Estimation::CellsDataContainer &container, size_t base_cell_ind) const
		{
			ul_list_t neighbour_cells = this->get_real_neighbour_cbs(container, base_cell_ind);
			if (neighbour_cells.empty())
				return -1;

			return this->get_best_merge_target(container, base_cell_ind, neighbour_cells);
		}

		long RealBarcodesMergeStrategy::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
															  const MergeStrategyAbstract::ul_list_t &neighbour_cells) const
		{
			if (neighbour_cells[0] == base_cell_ind)
				return base_cell_ind;

			double max_umigs_intersection_frac = 0;
			size_t best_neighbour_cell_ind = neighbour_cells[0];
			for (size_t neighbour_cell_ind: neighbour_cells)
			{
				size_t intersect_size = RealBarcodesMergeStrategy::get_umigs_intersect_size(container.cell_genes(base_cell_ind),
																							container.cell_genes(neighbour_cell_ind));
				double current_frac =  0.5 * intersect_size *( 1. / container.cell_size(base_cell_ind) + 1. / container.cell_size(neighbour_cell_ind));

				if (max_umigs_intersection_frac < current_frac)
				{
					max_umigs_intersection_frac = current_frac;
					best_neighbour_cell_ind = neighbour_cell_ind;
				}

//				container.stats().set(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
//									  container.cell_barcode(neighbour_cell_ind), current_frac);
			}

			if (max_umigs_intersection_frac < this->_min_merge_fraction)
				return -1;

			return best_neighbour_cell_ind;
		}

		RealBarcodesMergeStrategy::ul_list_t RealBarcodesMergeStrategy::get_real_neighbour_cbs(const Estimation::CellsDataContainer &container,
																							   size_t base_cell_ind) const
		{
			using BarcodesParsing::BarcodesParser;
			const std::string &base_cb = container.cell_barcode(base_cell_ind);
			std::vector<BarcodesParser::BarcodesDistance> barcodes_dists(this->_barcodes_parser->get_real_neighbour_cbs(base_cb));

			ul_list_t neighbour_cbs;
			if (barcodes_dists.empty())
				return neighbour_cbs;

			std::sort(barcodes_dists.begin(), barcodes_dists.end(),
			          [](const BarcodesParser::BarcodesDistance &d1, const BarcodesParser::BarcodesDistance &d2){
				          return d1.edit_distance < d2.edit_distance;});

			auto const &nearest_cb_parts = barcodes_dists.front();
			unsigned min_real_cb_dist = nearest_cb_parts.edit_distance;

			unsigned max_dist = this->get_max_merge_dist(min_real_cb_dist);
			for (auto const & cb_parts: barcodes_dists)
			{
				if (cb_parts.edit_distance > max_dist && !neighbour_cbs.empty())
					break;

				std::string cur_real_cb = this->_barcodes_parser->get_barcode(cb_parts.barcode_part_inds);
//				container.stats().set(Stats::MERGE_EDIT_DISTANCE_BY_CELL, base_cb, cur_real_cb, cb_parts.edit_distance);

				auto const current_cell_it = container.cell_ids_by_cb().find(cur_real_cb);
				if (current_cell_it != container.cell_ids_by_cb().end() &&
						container.cell_genes(current_cell_it->second).size() >= this->min_genes_before_merge() &&
						container.cell_size(current_cell_it->second) >= container.cell_size(base_cell_ind)) // Should pass equal sizes because it should pass current cell
				{
					neighbour_cbs.push_back(current_cell_it->second);
				}

				max_dist = std::max(max_dist, cb_parts.edit_distance); // If there are no cells on specified distance
			}

			return neighbour_cbs;
		}

		unsigned RealBarcodesMergeStrategy::get_max_merge_dist(unsigned min_real_cb_dist) const
		{
			return min_real_cb_dist;
		}

		std::string RealBarcodesMergeStrategy::merge_type() const
		{
			return "RealCBs";
		}
	}
}
