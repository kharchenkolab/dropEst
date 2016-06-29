#include "CellsDataContainer.h"

#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"
#include "Tools/UtilFunctions.h"

#include "api/BamMultiReader.h"

#include <boost/range/adaptor/reversed.hpp>
#include <limits>

using namespace std;

namespace Estimation
{
	const double EPS = 0.00001;

	CellsDataContainer::CellsDataContainer(double min_merge_fraction, int min_genes_before_merge, bool merge_tags,
	                                       int min_genes_after_merge, int max_merge_edit_distance, size_t top_print_size)
		: _merge_tags(merge_tags)
		, _top_print_size(top_print_size)
		, _min_merge_fraction(min_merge_fraction)
		, _min_genes_after_merge(min_genes_after_merge)
		, _min_genes_before_merge(min_genes_before_merge)
		, _max_merge_edit_distance(max_merge_edit_distance)
		, _is_initialized(false)
	{}

	void CellsDataContainer::merge_and_filter(const CellsDataContainer::s_ii_hash_t &umig_cells_counts)
	{
		if (!this->_is_initialized)
			throw runtime_error("You must initialize container");

		if (this->_merge_tags)
		{
			this->merge_cells(umig_cells_counts);

			this->filtered1_cells_genes_counts_sorted = this->count_filtered1_cells_genes(false);
		}
		else
		{
			for (auto const &gene_count : boost::adaptors::reverse(this->filtered1_cells_genes_counts_sorted))
			{
				if (gene_count.value < this->_min_genes_after_merge)
					break;

				this->_filtered_cells.push_back(gene_count.index);
			}
		}
	}

	void CellsDataContainer::merge_cells(const s_ii_hash_t &umig_cells_counts)
	{
		// cb merging
		int merges_count = 0;

		ISIHM cb_reassigned_to_it;
		ints_t cb_reassign_targets(this->_cells_genes.size());
		std::iota(cb_reassign_targets.begin(), cb_reassign_targets.end(), 0);

		L_TRACE << "merging linked tags ";

		int tag_index = 0;
		for (auto const &genes_count : this->filtered1_cells_genes_counts_sorted)
		{ // iterate through the minimally-selected CBs, from low to high counts
			if (++tag_index % 1000 == 0)
			{
				L_TRACE << "Total " << tag_index << " tags processed, " << merges_count << " cells merged";
			}

			i_i_hash_t umigs_intersect_top;
			size_t umigs_count = this->get_umigs_intersect_top(cb_reassign_targets, genes_count, umig_cells_counts,
			                                                   umigs_intersect_top);

			int top_cell_ind = -1;
			double top_cb_fraction = -1;
			long top_cb_genes_count = -1;
			for (auto const &cell: umigs_intersect_top)
			{
				int cell_ind = cell.first;
				double cb_fraction = cell.second / (double) umigs_count;
				if (cb_fraction - top_cb_fraction > EPS || (abs(cb_fraction - top_cb_fraction) < EPS &&
															  this->_cells_genes[cell_ind].size() > top_cb_genes_count))
				{
					top_cell_ind = cell_ind;
					top_cb_fraction = cb_fraction;
					top_cb_genes_count = this->_cells_genes[cell_ind].size();
				}
			}

			if (this->merge(top_cell_ind, top_cb_fraction, genes_count, cb_reassign_targets, cb_reassigned_to_it))
			{
				merges_count++;
			}
			else
			{
				this->_stats.add_merge_count(genes_count.value);
				if (genes_count.value >= this->_min_genes_after_merge)
				{
					this->_filtered_cells.push_back(genes_count.index);
				}
			}
		}

		this->_stats.merge(cb_reassign_targets, this->_cells_barcodes);

		if (this->_filtered_cells.size() > 1)
		{
			reverse(this->_filtered_cells.begin(), this->_filtered_cells.end());
		}

		L_INFO << "Done (" << merges_count << " merges performed)" << endl;
	}

	bool CellsDataContainer::merge(int target_cell_ind, double target_cell_fraction, const IndexedValue &source_genes_count,
	                               ints_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it)
	{
		if (target_cell_ind < 0)
			return false;

		// check if the top candidate is valid for merging
		if (target_cell_fraction < this->_min_merge_fraction)
			return false;

		size_t source_cell_ind = source_genes_count.index;
		int ed = Tools::edit_distance(this->_cells_barcodes[target_cell_ind].c_str(), this->_cells_barcodes[source_cell_ind].c_str());
		if (ed >= this->_max_merge_edit_distance)
			return false;

		merge_force(source_cell_ind, (size_t)target_cell_ind, source_genes_count.value, cb_reassign_targets, cb_reassigned_to_it);

		return true;
	}

	void CellsDataContainer::merge_force(size_t src_cell_id, size_t target_cell_ind, int target_genes_count,
	                                     CellsDataContainer::ints_t &cb_reassign_targets,
	                                     CellsDataContainer::ISIHM &cb_reassigned_to_it)
	{
		L_DEBUG << "Merge: " << _cells_barcodes[src_cell_id] << " to " << _cells_barcodes[target_cell_ind];
		this->stats().inc(Stats::MERGES_COUNT, _cells_barcodes[target_cell_ind]);

		_stats.add_merge_count(-1 * target_genes_count);

		for (auto const &gene: _cells_genes[src_cell_id])
		{
			for (auto const &umi_count: gene.second)
			{
				_cells_genes[target_cell_ind][gene.first][umi_count.first] += umi_count.second;
			}
		}

		this->reassign(src_cell_id, target_cell_ind, cb_reassign_targets, cb_reassigned_to_it);
	}

	void CellsDataContainer::reassign(size_t cell_id, size_t target_cell_id, ints_t &cb_reassign_targets,
	                                  ISIHM &cb_reassigned_to_it) const
	{
		cb_reassign_targets[cell_id] = target_cell_id; // set reassignment mapping
		cb_reassigned_to_it[target_cell_id].insert(cell_id); // reassign current cell

		// transfer mapping of the cbs previously mapped to kid
		auto reassigned_to_cell_iter = cb_reassigned_to_it.find(cell_id);
		if (reassigned_to_cell_iter == cb_reassigned_to_it.end())
			return;

		for (auto reassigned_id: reassigned_to_cell_iter->second)
		{
			cb_reassign_targets[reassigned_id] = target_cell_id; // update reassignment mapping
			cb_reassigned_to_it[target_cell_id].insert(reassigned_id);
		}

		reassigned_to_cell_iter->second.clear();
	}

	int CellsDataContainer::add_record(const string &cell_barcode, const string &umi, const string &gene, s_i_map_t &cells_ids)
	{
		if (this->_is_initialized)
			throw runtime_error("Container is already initialized");

		pair<s_i_map_t::iterator, bool> res = cells_ids.emplace(cell_barcode, this->_cells_genes.size());
		if (res.second)
		{ // new cb
			this->_cells_genes.push_back(genes_t());
			this->_cell_ids_by_cb[cell_barcode] = this->_cells_barcodes.size();
			this->_cells_barcodes.push_back(cell_barcode);
		}
		int cell_id = res.first->second;
		this->_cells_genes[cell_id][gene][umi]++;

		L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi] << " gene=" <<
				this->_cells_genes[cell_id][gene].size() << " CB=" << this->_cells_genes[cell_id].size();

		return cell_id;
	}

	CellsDataContainer::i_counter_t CellsDataContainer::count_filtered1_cells_genes(bool logs) const
	{
		i_counter_t cells_genes_counts; // <genes_count,cell_id> pairs
		for (size_t i = 0; i < this->_cells_genes.size(); i++)
		{
			if (this->_is_cell_excluded[i])
				continue;

			size_t genes_count = this->_cells_genes[i].size();
			if (genes_count >= this->_min_genes_before_merge)
			{
				cells_genes_counts.push_back(IndexedValue(i, genes_count));
			}
		}

		if (logs)
		{
			L_TRACE << cells_genes_counts.size() << " CBs with more than " << this->_min_genes_before_merge << " genes";
		}

		sort(cells_genes_counts.begin(), cells_genes_counts.end(), IndexedValue::value_less);

		if (logs)
		{
			L_TRACE << this->get_cb_count_top_verbose(cells_genes_counts);
		}

		return cells_genes_counts;
	}

	string CellsDataContainer::get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const
	{
		stringstream ss;
		if (cells_genes_counts.size() > 0)
		{
			ss << "top CBs:\n";
			size_t low_border = cells_genes_counts.size() - min(cells_genes_counts.size(), this->_top_print_size);
			for (size_t i = cells_genes_counts.size() - 1; i > low_border; --i)
			{
				ss << cells_genes_counts[i].value << "\t" << this->_cells_barcodes[cells_genes_counts[i].index] << "\n";
			}
		}
		else
		{
			ss << "no valid CBs found\n";
		}

		return ss.str();
	}

	Stats& CellsDataContainer::stats()
	{
		return this->_stats;
	}

	const Stats& CellsDataContainer::stats() const
	{
		return this->_stats;
	}

	const CellsDataContainer::genes_t &CellsDataContainer::cell_genes(size_t index) const
	{
		return this->_cells_genes[index];
	}

	const CellsDataContainer::i_counter_t &CellsDataContainer::cells_genes_counts_sorted() const
	{
		return this->filtered1_cells_genes_counts_sorted;
	}

	const string &CellsDataContainer::cell_barcode(size_t index) const
	{
		return this->_cells_barcodes[index];
	}

	const CellsDataContainer::ids_t &CellsDataContainer::filtered_cells() const
	{
		return this->_filtered_cells;
	}

	size_t CellsDataContainer::get_umigs_intersect_top(const ints_t &cb_reassign_targets,
	                                                   const IndexedValue &processed_genes_count,
	                                                   const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const
	{
		size_t umigs_count = 0;
		for (auto const &gene: this->_cells_genes[processed_genes_count.index])
		{
			const string &gene_name = gene.first;
			const s_i_map_t &umis = gene.second;

			for (auto const &umi_count: umis)
			{
				string umig = umi_count.first + gene_name;
				const i_i_hash_t &umig_cells = umigs_cells_counts.at(umig);
				for (auto const &cell : umig_cells)
				{
					int cell_with_same_umig_id = cell.first;
					long reassign_target = cb_reassign_targets[cell_with_same_umig_id]; //if not reassigned then cell_with_same_umig_id
					if (reassign_target == processed_genes_count.index)
						continue;

					if (this->_cells_genes[reassign_target].size() > processed_genes_count.value)
					{
						umig_top[reassign_target]++;
					}
				}
				umigs_count++;
			}
		}
		return umigs_count;
	}

	void CellsDataContainer::merge_by_real_barcodes(const string &barcodes_filename, size_t barcode2_length)
	{
		if (!this->_is_initialized)
			throw runtime_error("You must initialize container");

		names_t cbs1, cbs2;
		vector<bool> is_cell_real(this->_cells_genes.size(), false);

		CellsDataContainer::get_barcodes_list(barcodes_filename, cbs1, cbs2);
		if (cbs1.size() == 0)
			return;

		ISIHM cb_reassigned_to_it;
		ints_t cb_reassign_targets(this->_cells_genes.size());
		std::iota(cb_reassign_targets.begin(), cb_reassign_targets.end(), 0);

		size_t tag_index = 0, merges_count = 0;

		for (auto const &genes_count : boost::adaptors::reverse(this->filtered1_cells_genes_counts_sorted))
		{
			if (++tag_index % 1000 == 0)
			{
				L_TRACE << "Total " << tag_index << " tags processed, " << merges_count << " cells merged";
			}

			long real_cell_ind = this->get_real_cb(genes_count.index, cbs1, cbs2, barcode2_length);
			if (real_cell_ind == genes_count.index)
			{
				is_cell_real[genes_count.index] = true;
				continue;
			}

			if (real_cell_ind < 0)
			{
				this->_is_cell_excluded[genes_count.index] = true;
				continue;
			}

			this->merge_force(genes_count.index, real_cell_ind, genes_count.value, cb_reassign_targets, cb_reassigned_to_it);
			merges_count++;
		}
		L_INFO << "Total " << merges_count << " merges";

		this->filtered1_cells_genes_counts_sorted = this->count_filtered1_cells_genes(true);
		for (auto const &gene_count : boost::adaptors::reverse(this->filtered1_cells_genes_counts_sorted))
		{
			if (gene_count.value < this->_min_genes_after_merge)
				break;

			if (!is_cell_real[gene_count.index])
				continue;

			L_DEBUG << "Add cell to filtered: " << gene_count.value << " " << gene_count.index;
			this->_filtered_cells.push_back(gene_count.index);
		}

		this->_stats.merge(cb_reassign_targets, this->_cells_barcodes);
	}

	long CellsDataContainer::get_real_cb(size_t base_cell_ind, const names_t &cbs1, const names_t &cbs2,
	                                            size_t barcode2_length) const
	{
		string base_cb = this->_cells_barcodes[base_cell_ind];
		i_counter_t dists1, dists2;

		string cb_part1 = base_cb.substr(0, base_cb.length() - barcode2_length);
		string cb_part2 = base_cb.substr(base_cb.length() - barcode2_length + 1);

		CellsDataContainer::fill_distances_to_cb(cb_part1, cb_part2, cbs1, cbs2, dists1, dists2);

		if (dists1[0].value == 0 && dists2[0].value == 0)
			return base_cell_ind;

		L_DEBUG <<"Get real neighbours to " << base_cb;
		ids_t neighbour_cells = this->get_real_neighbour_cbs(cbs1, cbs2, base_cb, dists1, dists2);
		if (neighbour_cells.empty())
			return -1;

		size_t max_umigs_intersection_size = 0;
		size_t best_neighbour_cell_ind = neighbour_cells[0];
		for (size_t neighbour_cell_ind: neighbour_cells)
		{
			size_t umigs_intersection_size = this->get_umigs_intersection_size(base_cell_ind, neighbour_cell_ind);
			if (max_umigs_intersection_size < umigs_intersection_size)
			{
				max_umigs_intersection_size = umigs_intersection_size;
				best_neighbour_cell_ind = neighbour_cell_ind;
			}
			this->_stats.add_str(Stats::MERGE_INTERSECT_SIZE_BY_CELL, this->_cells_barcodes[neighbour_cell_ind],
			                     base_cb, umigs_intersection_size);
		}
		this->_stats.add_str(Stats::MERGE_REAL_INTERSECT_SIZE_BY_CELL, this->_cells_barcodes[best_neighbour_cell_ind],
		                     base_cb, max_umigs_intersection_size);

		return best_neighbour_cell_ind;
	}

	void CellsDataContainer::fill_distances_to_cb(const string &cb_part1, const string &cb_part2, const names_t &cbs1,
	                                              const names_t &cbs2, CellsDataContainer::i_counter_t &dists1,
	                                              CellsDataContainer::i_counter_t &dists2)
	{
		dists1.reserve(cbs1.size());
		dists2.reserve(cbs2.size());

		size_t i = 0;
		for (auto const cb1: cbs1)
		{
			dists1.push_back(CellsDataContainer::IndexedValue(i, Tools::edit_distance(cb_part1.c_str(), cb1.c_str())));
			i++;
		}

		i = 0;
		for (auto const cb2: cbs2)
		{
			dists2.push_back(CellsDataContainer::IndexedValue(i, Tools::edit_distance(cb_part2.c_str(), cb2.c_str())));
			i++;
		}

		sort(dists1.begin(), dists1.end(), IndexedValue::value_less);
		sort(dists2.begin(), dists2.end(), IndexedValue::value_less);
	}

	/*
	 * Return the list of cbs with minimal distances that were found in the current dataset (work not really precisely)
	 */
	/*CellsDataContainer::ids_t CellsDataContainer::get_real_neighbour_cbs(const names_t &cbs1, const names_t &cbs2,
	                                                   const string &base_cb, const CellsDataContainer::i_counter_t &dists1,
	                                                   const CellsDataContainer::i_counter_t &dists2) const
	{
		ids_t neighbour_cbs;
		long prev_dist1 = std::numeric_limits<long>::max();
		for (size_t ind1 = 0; ind1 < dists1.size(); ++ind1)
		{
			const IndexedValue &cur_dist1 = dists1[ind1];
			if (cur_dist1.value > prev_dist1 && !neighbour_cbs.empty())
				break;

			if (neighbour_cbs.empty() && cur_dist1.value > prev_dist1)
			{
				L_DEBUG << "There is no neighbours, cb1, prev: " << prev_dist1 << ", cur: " << cur_dist1.value;
			}

			long prev_dist2 = std::numeric_limits<long>::max();
			for (size_t ind2 = 0; ind2 < dists2.size(); ++ind2)
			{
				const IndexedValue &cur_dist2 = dists2[ind2];
				if (cur_dist2.value > prev_dist2 &&
						(!neighbour_cbs.empty() ||
								ind1 != dists1.size() - 1 && ind2 != dists2.size() - 1 &&
								dists2[ind2 + 1].value - dists2[0].value > dists1[ind2 + 1].value - dists1[0].value))
					break;

				if (neighbour_cbs.empty() && cur_dist2.value > prev_dist2)
				{
					L_DEBUG << "There is no neighbours, cb2, prev: " << prev_dist2 << ", cur: " << cur_dist2.value;
				}

				string current_cb = cbs1[cur_dist1.index] + cbs2[cur_dist2.index];
				auto const current_cell_it = this->_cell_ids_by_cb.find(current_cb);
				if (current_cell_it != this->_cell_ids_by_cb.end())
				{
					neighbour_cbs.push_back(current_cell_it->second);
					this->_stats.add_str(Stats::MERGE_EDIT_DISTANCE_BY_CELL, current_cb, base_cb, cur_dist1.value + cur_dist2.value);
				}
				else
				{
					this->_stats.add_str(Stats::MERGE_REJECTION_BY_CELL, current_cb, base_cb, cur_dist1.value + cur_dist2.value);
				}

				prev_dist2 = cur_dist2.value;
			}

			if (!neighbour_cbs.empty() && (ind1 == dists1.size() || dists1[ind1 + 1].value >= cur_dist1.value))
				break;

			prev_dist1 = cur_dist1.value;
		}

		return neighbour_cbs;
	}*/

	CellsDataContainer::ids_t CellsDataContainer::get_real_neighbour_cbs(const names_t &cbs1, const names_t &cbs2,
	                                                                     const string &base_cb, const CellsDataContainer::i_counter_t &dists1,
	                                                                     const CellsDataContainer::i_counter_t &dists2) const
	{
		typedef tuple<size_t, size_t, long> tuple_t;
		vector<tuple_t> real_cell_inds;
		for (size_t ind1 = 0; ind1 < dists1.size(); ++ind1)
		{
			auto const &cur_dist1 = dists1[ind1];
			if (cur_dist1.value + dists2[0].value > CellsDataContainer::MAX_REAL_MERGE_EDIT_DISTANCE)
				break;
			
			for (size_t ind2 = 0; ind2 < dists2.size(); ++ind2)
			{
				auto const &cur_dist2 = dists2[ind2];
				long cur_ed = cur_dist1.value + cur_dist2.value;
				if (cur_ed > CellsDataContainer::MAX_REAL_MERGE_EDIT_DISTANCE)
					break;

				real_cell_inds.push_back(make_tuple(cur_dist1.index, cur_dist2.index, cur_ed));
			}
		}

		sort(real_cell_inds.begin(), real_cell_inds.end(), [](const tuple_t &t1, const tuple_t &t2){ return get<2>(t1) < get<2>(t2);});

		ids_t neighbour_cbs;
		long prev_dist = std::numeric_limits<long>::max();
		for (auto const & inds: real_cell_inds)
		{
			long cur_ed = get<2>(inds);
			if (cur_ed > prev_dist && !neighbour_cbs.empty())
				break;

			string current_cb = cbs1[get<0>(inds)] + cbs2[get<1>(inds)];
			auto const current_cell_it = this->_cell_ids_by_cb.find(current_cb);
			if (current_cell_it != this->_cell_ids_by_cb.end())
			{
				neighbour_cbs.push_back(current_cell_it->second);
				this->_stats.add_str(Stats::MERGE_EDIT_DISTANCE_BY_CELL, current_cb, base_cb, cur_ed);
			}
			else
			{
				this->_stats.add_str(Stats::MERGE_REJECTION_BY_CELL, current_cb, base_cb, cur_ed);
			}
			prev_dist = cur_ed;
		}

		return neighbour_cbs;
	}

	void CellsDataContainer::get_barcodes_list(const std::string &barcodes_filename,
	                                           std::vector<std::string> &barcodes1,
	                                           std::vector<std::string> &barcodes2)
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

		if (barcodes1.size() == 0)
		{
			L_WARN << "WARNING: empty barcodes list";
		}
	}

	size_t CellsDataContainer::get_umigs_intersection_size(size_t cell1_ind, size_t cell2_ind) const
	{
		const auto &cell1 = this->_cells_genes[cell1_ind];
		const auto &cell2 = this->_cells_genes[cell2_ind];

		auto gene1_it = cell1.begin();
		auto gene2_it = cell2.begin();

		size_t intersect_size = 0;
		while (gene1_it != cell1.end() && gene2_it != cell2.end())
		{
			int comp_res = gene1_it->first.compare(gene2_it->first);
			if (comp_res < 0)
			{
				gene1_it++;
				continue;
			}

			if (comp_res > 0)
			{
				gene2_it++;
				continue;
			}

			auto umi1_it = gene1_it->second.begin();
			auto umi2_it = gene2_it->second.begin();

			while (umi1_it != gene1_it->second.end() && umi2_it != gene2_it->second.end())
			{
				comp_res = umi1_it->first.compare(umi2_it->first);
				if (comp_res < 0)
				{
					umi1_it++;
					continue;
				}

				if (comp_res > 0)
				{
					umi2_it++;
					continue;
				}

				++intersect_size;
				umi2_it++;
				umi1_it++;
			}

			gene1_it++;
			gene2_it++;
		}

		return intersect_size;
	}

	void CellsDataContainer::set_initialized()
	{
		this->_is_cell_excluded = flags_t(this->_cells_barcodes.size(), false);
		this->filtered1_cells_genes_counts_sorted = this->count_filtered1_cells_genes();
		this->_is_initialized = true;
	}

	CellsDataContainer::names_t CellsDataContainer::excluded_cells() const
	{
		names_t ec;
		for (size_t i = 0; i < this->_is_cell_excluded.size(); ++i)
		{
			if (this->_is_cell_excluded[i])
			{
				ec.push_back(this->_cells_barcodes[i]);
			}
		}

		return ec;
	}
}
