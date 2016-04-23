#pragma once

#include "Stats.h"

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <vector>

namespace TestEstimator
{
	struct testMerge;
}

namespace Estimation
{
	class CellsDataContainer
	{
		friend struct TestEstimator::testMerge;

	private:
		typedef boost::unordered_set<size_t> i_set_t;
		typedef boost::unordered_map<int, int> i_i_hash_t;
		typedef boost::unordered_map<size_t, i_set_t> ISIHM;

	public:
		struct IndexedCount
		{
			IndexedCount(size_t cell_index, size_t count)
					: cell_index(cell_index), count(count)
			{}

			size_t cell_index;
			size_t count;

			static bool counts_comp(const IndexedCount &ic1, const IndexedCount &ic2)
			{
				return ic1.count < ic2.count;
			}
		};

		typedef boost::unordered_map<std::string, int> s_i_hash_t;
		typedef boost::unordered_map<std::string, i_i_hash_t> s_ii_hash_t;
		typedef boost::unordered_map<std::string, s_i_hash_t> genes_t;

		typedef std::vector<IndexedCount> i_counter_t;
		typedef std::vector<long> ints_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<std::string> names_t;

	private:
		const bool _merge_tags;
		const size_t _top_print_size;
		const double _min_merge_fraction;
		const int _min_genes_before_merge;
		const int _min_genes_after_merge;
		const int _max_merge_edit_distance;

		std::vector<genes_t> _cells_genes; //cell_id -> gen_name -> umi -> count
		names_t _cells_barcodes;
		ids_t _filtered_cells;

		i_counter_t _cells_genes_counts_sorted;

		Stats _stats;

	private:
		void merge_genes(const s_ii_hash_t &umig_cbs);

		bool merge(int top_cell_ind, double top_cell_fraction, const IndexedCount &processed_gene_count,
				   ints_t &cb_reassigned, ISIHM &cb_reassigned_to);

		void reassign(size_t cell_id, int target_cell_id, ints_t &cb_reassigned, ISIHM &cb_reassigned_to) const;

		size_t get_umig_top(const ints_t &cb_reassigned, const IndexedCount &processed_genes_count,
							const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const;

		i_counter_t count_cells_genes(bool logs = true) const;

		std::string get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const;

	public:
		CellsDataContainer(double min_merge_fraction, int min_genes_before_merge, bool merge_tags,
						   int min_genes_after_merge, int max_merge_edit_distance, size_t top_print_size);

		void init_filled(const CellsDataContainer::s_ii_hash_t &umig_cells_counts);

		int add_record(const std::string &cell_barcode, const std::string &umi, const std::string &gene, s_i_hash_t &cells_ids);

		Stats &stats();
		const Stats &stats() const;

		const ids_t &filtered_cells() const;

		const i_counter_t &cells_genes_counts_sorted() const;

		const genes_t &cell_genes(size_t index) const;

		const std::string &cell_barcode(size_t index) const;
	};
}