#pragma once

#include "Stats.h"

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <map>
#include <vector>

namespace TestEstimator
{
	struct testMerge;
	struct testBarcodesFile;
	struct testUmigsIntersection;
	struct testFillDistances;
	struct testRealNeighboursCbs;
	struct testRealNeighbours;
}

namespace Estimation
{
	class CellsDataContainer
	{
		friend struct TestEstimator::testMerge;
		friend struct TestEstimator::testBarcodesFile;
		friend struct TestEstimator::testUmigsIntersection;
		friend struct TestEstimator::testFillDistances;
		friend struct TestEstimator::testRealNeighboursCbs;
		friend struct TestEstimator::testRealNeighbours;

	private:
		typedef boost::unordered_set<size_t> i_set_t;
		typedef boost::unordered_map<int, int> i_i_hash_t;
		typedef boost::unordered_map<size_t, i_set_t> ISIHM;

	public:
		struct IndexedValue
		{
			IndexedValue(size_t cell_index, size_t count)
					: index(cell_index), value(count)
			{}

			size_t index;
			long value;

			static bool value_less(const IndexedValue &ic1, const IndexedValue &ic2)
			{
				return ic1.value < ic2.value;
			}
		};

		typedef boost::unordered_map<std::string, size_t> s_ul_hash_t;
		typedef boost::unordered_map<std::string, i_i_hash_t> s_ii_hash_t;

		typedef std::map<std::string, int> s_i_map_t;
		typedef std::map<std::string, s_i_map_t> genes_t;

		typedef std::vector<IndexedValue> i_counter_t;
		typedef std::vector<long> ints_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<std::string> names_t;
		typedef std::vector<bool> flags_t;

	private:
		static const int MAX_REAL_MERGE_EDIT_DISTANCE=5;
		const bool _merge_tags;
		const size_t _top_print_size;
		const double _min_merge_fraction;
		const int _min_genes_before_merge;
		const int _min_genes_after_merge;
		const int _max_merge_edit_distance;

		std::vector<genes_t> _cells_genes; //cell_id -> gen_name -> umi -> count
		names_t _cells_barcodes;
		flags_t _is_cell_excluded;
		s_ul_hash_t _cell_ids_by_cb;
		ids_t _filtered_cells;

		i_counter_t filtered1_cells_genes_counts_sorted;
		bool _is_initialized;

		mutable Stats _stats;

	private:
		void merge_cells(const s_ii_hash_t &umig_cbs);

		bool merge(int target_cell_ind, double target_cell_fraction, const IndexedValue &source_genes_count,
				   ints_t &cb_reassigned, ISIHM &cb_reassigned_to);

		void reassign(size_t cell_id, size_t target_cell_id, ints_t &cb_reassigned, ISIHM &cb_reassigned_to) const;

		size_t get_umigs_intersect_top(const ints_t &cb_reassigned, const IndexedValue &processed_genes_count,
		                               const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const;

		i_counter_t count_filtered1_cells_genes(bool logs = true) const;

		ids_t get_real_neighbour_cbs(const names_t &cbs1, const names_t &cbs2, const std::string &base_cb,
		                             const i_counter_t &dists1, const i_counter_t &dists2) const;

		long get_real_cb(size_t base_cell_ind, const names_t &cbs1, const names_t &cbs2, size_t barcode2_length) const;

		size_t get_umigs_intersection_size(size_t cell1_ind, size_t cell2_ind) const;

		void merge_force(size_t src_cell_id, size_t target_cell_ind, int processed_gene_count,
		                 ints_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it);

		std::string get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const;

		static void get_barcodes_list(const std::string &barcodes_filename, std::vector<std::string> &barcodes1,
		                              std::vector<std::string> &barcodes2);

		static void fill_distances_to_cb(const std::string &cb_part1, const std::string &cb_part2, const names_t &cbs1,
		                                 const names_t &cbs2, i_counter_t &dists1, i_counter_t &dists2);

	public:
		CellsDataContainer(double min_merge_fraction, int min_genes_before_merge, bool merge_tags,
						   int min_genes_after_merge, int max_merge_edit_distance, size_t top_print_size);

		void merge_and_filter(const CellsDataContainer::s_ii_hash_t &umig_cells_counts);

		void merge_by_real_barcodes(const std::string &barcodes_filename, size_t barcode2_length);

		int add_record(const std::string &cell_barcode, const std::string &umi, const std::string &gene, s_i_map_t &cells_ids);

		Stats &stats();
		const Stats &stats() const;

		const ids_t &filtered_cells() const;

		const i_counter_t &cells_genes_counts_sorted() const;

		const genes_t &cell_genes(size_t index) const;

		const std::string &cell_barcode(size_t index) const;

		names_t excluded_cells() const;

		void set_initialized();
	};
}