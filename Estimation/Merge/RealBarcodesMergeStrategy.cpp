#include <boost/range/adaptor/reversed.hpp>
#include <Tools/Logs.h>
#include <fstream>
#include <Tools/UtilFunctions.h>
#include "RealBarcodesMergeStrategy.h"

#include <omp.h>

namespace Estimation
{
	namespace Merge
	{
		using Tools::IndexedValue;

		RealBarcodesMergeStrategy::RealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
															 int min_genes_before_merge, int min_genes_after_merge,
															 int max_merge_edit_distance, double min_merge_fraction)
				: MergeStrategyBase(min_genes_before_merge, min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
				, _barcodes_filename(barcodes_filename)
				, _barcode2_length(barcode2_length)
		{}

		void RealBarcodesMergeStrategy::get_barcodes_list(const std::string &barcodes_filename, names_t &barcodes1,
														  names_t &barcodes2)
		{
			std::ifstream cb_f(barcodes_filename);
			if (cb_f.fail())
				throw std::runtime_error("Can't open barcodes file: '" + barcodes_filename + "'");

			std::string line;
			while (std::getline(cb_f, line))
			{
				size_t space_ind = line.find(' ');
				if (space_ind == std::string::npos)
				{
					L_WARN << "WARNING: barcodes line has bad format: '" << line << "'";
					continue;
				}

				barcodes1.push_back(Tools::reverse_complement(line.substr(0, space_ind)));
				barcodes2.push_back(Tools::reverse_complement(line.substr(space_ind + 1)));
			}
		}

		long RealBarcodesMergeStrategy::get_merge_target(const Estimation::CellsDataContainer &container, size_t base_cell_ind) const
		{
			std::string base_cb = container.cell_barcode(base_cell_ind);
			i_counter_t dists1, dists2;

			std::string cb_part1 = base_cb.substr(0, base_cb.length() - this->_barcode2_length);
			std::string cb_part2 = base_cb.substr(base_cb.length() - this->_barcode2_length);

			this->fill_distances_to_cb(cb_part1, cb_part2, dists1, dists2);

//			L_DEBUG <<"Get real neighbours to " << base_cb;
			ul_list_t neighbour_cells = this->get_real_neighbour_cbs(container, base_cell_ind, dists1, dists2);
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

				#pragma omp critical(MERGE_PROB_BY_CELL)
				container.stats().set(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
									  container.cell_barcode(neighbour_cell_ind), current_frac);
			}

			if (max_umigs_intersection_frac < this->_min_merge_fraction)
				return -1;

			return best_neighbour_cell_ind;
		}

		RealBarcodesMergeStrategy::ul_list_t RealBarcodesMergeStrategy::get_real_neighbour_cbs(const Estimation::CellsDataContainer &container,
																							   size_t base_cell_ind,
																							   const i_counter_t &dists1,
																							   const i_counter_t &dists2) const
		{
			const std::string &base_cb = container.cell_barcode(base_cell_ind);
			typedef std::tuple<size_t, size_t, long> tuple_t;
			std::vector<tuple_t> real_cell_inds;
			for (size_t ind1 = 0; ind1 < dists1.size(); ++ind1)
			{
				auto const &cur_dist1 = dists1[ind1];
				if (cur_dist1.value + dists2[0].value > RealBarcodesMergeStrategy::MAX_REAL_MERGE_EDIT_DISTANCE)
					break;

				for (size_t ind2 = 0; ind2 < dists2.size(); ++ind2)
				{
					auto const &cur_dist2 = dists2[ind2];
					long cur_ed = cur_dist1.value + cur_dist2.value;
					if (cur_ed > RealBarcodesMergeStrategy::MAX_REAL_MERGE_EDIT_DISTANCE)
						break;

					real_cell_inds.push_back(std::make_tuple(cur_dist1.index, cur_dist2.index, cur_ed));
				}
			}

			ul_list_t neighbour_cbs;
			if (real_cell_inds.empty())
				return neighbour_cbs;

			sort(real_cell_inds.begin(), real_cell_inds.end(), [](const tuple_t &t1, const tuple_t &t2){ return std::get<2>(t1) < std::get<2>(t2);});

			long max_dist = this->get_max_merge_dist(std::get<2>(real_cell_inds.front()));
			for (auto const & inds: real_cell_inds)
			{
				long cur_ed = std::get<2>(inds);
				if (cur_ed > max_dist && !neighbour_cbs.empty())
					break;

				std::string current_cb = this->_barcodes1[std::get<0>(inds)] + this->_barcodes2[std::get<1>(inds)];
				auto const current_cell_it = container.cell_ids_by_cb().find(current_cb);
				if (current_cell_it != container.cell_ids_by_cb().end() &&
						container.cell_genes(current_cell_it->second).size() >= this->min_genes_before_merge() &&
						container.cell_size(current_cell_it->second) > container.cell_size(base_cell_ind))
				{
					neighbour_cbs.push_back(current_cell_it->second);
					container.stats().add(Estimation::Stats::MERGE_EDIT_DISTANCE_BY_CELL, current_cb, base_cb, cur_ed);
				}
				else
				{
					container.stats().add(Estimation::Stats::MERGE_REJECTION_BY_CELL, current_cb, base_cb, cur_ed);
				}

				max_dist = std::max(max_dist, cur_ed); // If there are no cells on specified distance
			}

			return neighbour_cbs;
		}

		void RealBarcodesMergeStrategy::fill_distances_to_cb(const std::string &cb_part1, const std::string &cb_part2,
															 i_counter_t &dists1, i_counter_t &dists2) const
		{
			dists1.reserve(this->_barcodes1.size());
			dists2.reserve(this->_barcodes2.size());

			size_t i = 0;
			for (auto const cb1: this->_barcodes1)
			{
				dists1.push_back(IndexedValue(i, Tools::edit_distance(cb_part1.c_str(), cb1.c_str())));
				i++;
			}

			i = 0;
			for (auto const cb2: this->_barcodes2)
			{
				dists2.push_back(IndexedValue(i, Tools::edit_distance(cb_part2.c_str(), cb2.c_str())));
				i++;
			}

			sort(dists1.begin(), dists1.end(), IndexedValue::value_less);
			sort(dists2.begin(), dists2.end(), IndexedValue::value_less);
		}

		void RealBarcodesMergeStrategy::init(const Estimation::CellsDataContainer &container)
		{
			MergeStrategyAbstract::init(container);
			RealBarcodesMergeStrategy::get_barcodes_list(this->_barcodes_filename, this->_barcodes1, this->_barcodes2);

			if (this->_barcodes1.size() == 0)
				throw std::runtime_error("ERROR: empty barcodes list");
		}

		void RealBarcodesMergeStrategy::release()
		{
			this->_barcodes1.clear();
			this->_barcodes2.clear();
			MergeStrategyAbstract::release();
		}

		long RealBarcodesMergeStrategy::get_max_merge_dist(long min_real_cb_dist) const
		{
			return min_real_cb_dist;
		}

		bool RealBarcodesMergeStrategy::need_parallel() const
		{
			return false;
		}

		std::string RealBarcodesMergeStrategy::merge_type() const
		{
			return "RealCBs";
		}
	}
}
