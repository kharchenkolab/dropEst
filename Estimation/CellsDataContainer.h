#pragma once

#include "Stats.h"
#include "Tools/IndexedValue.h"

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <map>
#include <vector>

namespace TestEstimator
{
	struct testMerge;
}

namespace Estimation
{
	namespace Merge
	{
		class IMergeStrategy;
	}

	class CellsDataContainer
	{
		friend struct TestEstimator::testMerge;

	public:
		typedef boost::unordered_map<int, int> i_i_hash_t;
		typedef boost::unordered_map<std::string, size_t> s_ul_hash_t;
		typedef boost::unordered_map<std::string, i_i_hash_t> s_ii_hash_t;

		typedef std::map<std::string, int> s_i_map_t;
		typedef std::map<std::string, s_i_map_t> genes_t;

		typedef std::vector<Tools::IndexedValue> i_counter_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<std::string> names_t;
		typedef std::vector<bool> flags_t;

	private:
		std::shared_ptr<Merge::IMergeStrategy> _merge_strategy;

		const size_t _top_print_size;

		std::vector<genes_t> _cells_genes; //cell_id -> gen_name -> umi -> count
		names_t _cell_barcodes;
		flags_t _is_cell_excluded;
		s_ul_hash_t _cell_ids_by_cb;
		ids_t _filtered_cells;

		i_counter_t filtered_cells_genes_counts_sorted;
		bool _is_initialized;

		mutable Stats _stats;

	private:
		std::string get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const;

	public:
		CellsDataContainer(std::shared_ptr<Merge::IMergeStrategy> merge_strategy, size_t top_print_size);

		void merge_and_filter(const CellsDataContainer::s_ii_hash_t &umig_cells_counts);

		int add_record(const std::string &cell_barcode, const std::string &umi, const std::string &gene, s_i_map_t &cells_ids);

		void merge(size_t source_cell_ind, size_t target_cell_ind);

		void exclude_cell(size_t index);

		void update_cells_genes_counts(int threshold, bool logs = true);

		const s_ul_hash_t& cell_ids_by_cb() const;

		Stats &stats() const;

		const ids_t &filtered_cells() const;

		const i_counter_t &cells_genes_counts_sorted() const;

		const genes_t &cell_genes(size_t index) const;

		const std::string &cell_barcode(size_t index) const;

		names_t excluded_cells() const;

		void set_initialized();

		const names_t &cell_barcodes() const;
	};
}