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
		class MergeStrategyAbstract;
	}

	class CellsDataContainer
	{
		friend struct TestEstimator::testMerge;

	public:
		typedef boost::unordered_map<std::string, size_t> s_ul_hash_t;

		typedef std::map<std::string, size_t> s_i_map_t;
		typedef std::map<std::string, s_i_map_t> genes_t;

		typedef std::vector<Tools::IndexedValue> i_counter_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<size_t> counts_t;
		typedef std::vector<std::string> names_t;
		typedef std::vector<bool> flags_t;

	private:
		std::shared_ptr<Merge::MergeStrategyAbstract> _merge_strategy;

		const size_t _top_print_size;

		std::vector<genes_t> _cells_genes; //cell_id -> gen_name -> umi -> count
		names_t _cell_barcodes;
		flags_t _is_cell_excluded;
		flags_t _is_cell_merged;
		s_ul_hash_t _cell_ids_by_cb;
		ids_t _filtered_cells;
		ids_t _merge_targets;

		counts_t _cell_sizes;
		i_counter_t _filtered_cells_gene_counts_sorted;
		bool _is_initialized;

		mutable Stats _stats;

	private:
		std::string get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const;

	public:
		CellsDataContainer(std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy, size_t top_print_size);

		void merge_and_filter();

		size_t add_record(const std::string &cell_barcode, const std::string &umi, const std::string &gene);

		void merge(size_t source_cell_ind, size_t target_cell_ind);

		void exclude_cell(size_t index);

		void update_cell_sizes(int genes_threshold, bool logs = true);

		const s_ul_hash_t& cell_ids_by_cb() const;

		Stats &stats() const;

		const ids_t &filtered_cells() const;

		const i_counter_t &cells_gene_counts_sorted() const;

		const genes_t &cell_genes(size_t index) const;

		const std::string &cell_barcode(size_t index) const;

		names_t excluded_cells() const;

		bool is_cell_merged(size_t cell_id) const;
		bool is_cell_excluded(size_t cell_id) const;

		const ids_t& merge_targets() const;

		void set_initialized();

		const names_t &cell_barcodes_raw() const;

		size_t const cell_size(size_t cell_index) const;

		s_ul_hash_t umis_distribution() const;

		std::string merge_type() const;
	};
}