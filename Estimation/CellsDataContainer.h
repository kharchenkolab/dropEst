#pragma once

#include "Stats.h"
#include "Tools/IndexedValue.h"

#include <string>

#include <map>
#include <vector>
#include <unordered_map>

namespace TestEstimator
{
	struct testMerge;
	struct testGeneMatchLevelUmiExclusion;
}

namespace Estimation
{
	namespace Merge
	{
		class MergeStrategyAbstract;
	}

	namespace MergeUMIs
	{
		class MergeUMIsStrategySimple;
	}

	class CellsDataContainer
	{
		friend struct TestEstimator::testMerge;
		friend struct TestEstimator::testGeneMatchLevelUmiExclusion;

	public:
		typedef std::unordered_map<std::string, size_t> s_ul_hash_t;
		typedef std::unordered_map<std::string, std::string> s_s_hash_t;

		typedef long umi_cnt_t;
		typedef std::map<std::string, umi_cnt_t> s_i_map_t;
		typedef std::map<std::string, s_i_map_t> genes_t;

		typedef std::vector<Tools::IndexedValue> i_counter_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<size_t> counts_t;
		typedef std::vector<std::string> names_t;
		typedef std::vector<bool> flags_t;

	private:
		std::shared_ptr<Merge::MergeStrategyAbstract> _merge_strategy;
		std::shared_ptr<MergeUMIs::MergeUMIsStrategySimple> _umi_merge_strategy;

		const size_t _top_print_size;
		static const umi_cnt_t UMI_EXCLUDED = -100000;

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
		void remove_excluded_umis();

		std::string get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const;

	public:
		CellsDataContainer(std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy,
		                   std::shared_ptr<MergeUMIs::MergeUMIsStrategySimple> umi_merge_strategy, size_t top_print_size);

		size_t add_record(const std::string &cell_barcode, const std::string &umi, const std::string &gene, bool set_excluded=false);
		void exclude_cell(size_t index);

		void merge_and_filter();
		void merge_cells(size_t source_cell_ind, size_t target_cell_ind);

		void merge_umis(size_t cell_id, const std::string &gene, const s_s_hash_t &merge_targets);

		void update_cell_sizes(int genes_threshold, bool logs = true);

		void set_initialized();

		const names_t &cell_barcodes_raw() const;
		const i_counter_t &cells_gene_counts_sorted() const;
		const s_ul_hash_t& cell_ids_by_cb() const;
		names_t excluded_cells() const;
		const ids_t &filtered_cells() const;
		const ids_t& merge_targets() const;

		const std::string &cell_barcode(size_t index) const;
		const genes_t &cell_genes(size_t index) const;
		size_t cell_size(size_t cell_index) const;
		bool is_cell_excluded(size_t cell_id) const;
		bool is_cell_merged(size_t cell_id) const;

		Stats &stats() const;
		std::string merge_type() const;

		s_ul_hash_t umis_distribution() const;
	};
}