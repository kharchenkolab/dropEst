#pragma once

#include "Cell.h"
#include "Stats.h"
#include "Tools/IndexedValue.h"
#include <Tools/GtfRecord.h>
#include "UMI.h"
#include "StringIndexer.h"

#include <string>

#include <map>
#include <vector>
#include <unordered_map>
#include <fstream>

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
		namespace UMIs
		{
			class MergeUMIsStrategySimple;
		}
	}

	class CellsDataContainer
	{
		friend struct TestEstimator::testMerge;
		friend struct TestEstimator::testGeneMatchLevelUmiExclusion;

	public:
		typedef std::unordered_map<std::string, std::string> s_s_hash_t;
		typedef std::unordered_map<std::string, size_t> s_ul_hash_t;
		typedef std::unordered_map<std::string, int> s_i_hash_t; // not long because of RCpp

		typedef std::vector<Tools::IndexedValue> i_counter_t;
		typedef std::vector<size_t> ids_t;
		typedef std::vector<int> counts_t;
		typedef std::vector<std::string> names_t;

	private:
		std::shared_ptr<Merge::MergeStrategyAbstract> _merge_strategy;
		std::shared_ptr<Merge::UMIs::MergeUMIsStrategySimple> _umi_merge_strategy;

		static const size_t TOP_PRINT_SIZE;
		const int _max_cells_num;

		std::vector<Cell> _cells; //cell_id -> gen_name -> umi -> #reads
		s_ul_hash_t _cell_ids_by_cb;
		ids_t _filtered_cells; // Sorted ascending by the number of genes
		ids_t _merge_targets;

		bool _is_initialized;
		const UMI::Mark::query_t _query_marks;

		size_t _has_exon_reads;
		size_t _has_intron_reads;
		size_t _has_not_annotated_reads;
		size_t _number_of_real_cells;

		StringIndexer _umi_indexer;
		StringIndexer _gene_indexer;

	private:
		std::string get_cb_count_top_verbose() const;
		size_t update_cell_sizes(const UMI::Mark::query_t &query_marks, size_t requested_genes_threshold, int cell_threshold);
		void update_cell_stats(size_t cell_id, const UMI::Mark &mark, const std::string &chromosome_name);

		bool compare_cells(size_t cell1_id, size_t cell2_id) const;

		size_t update_filtered_gene_counts(size_t requested_genes_threshold, int cell_threshold);

	public:
		CellsDataContainer(std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy,
		                   std::shared_ptr<Merge::UMIs::MergeUMIsStrategySimple> umi_merge_strategy,
		                   const std::vector<UMI::Mark> &gene_match_levels, int max_cells_num = -1);

		void add_record(const std::string &cell_barcode, const std::string &umi, const std::string &gene,
		                const std::string &chr_name = "", const UMI::Mark &umi_mark = UMI::Mark::HAS_EXONS);
		void exclude_cell(size_t index);

		void merge_and_filter();
		void merge_cells(size_t source_cell_ind, size_t target_cell_ind);

		void merge_umis(size_t cell_id, StringIndexer::index_t gene, const s_s_hash_t &merge_targets);

		void set_initialized();

		size_t total_cells_number() const;
		size_t cell_id_by_cb(const std::string &barcode) const;
		const ids_t& filtered_cells() const;
		const ids_t& merge_targets() const;
		const UMI::Mark::query_t& gene_match_level() const;

		s_i_hash_t get_stat_by_real_cells(Stats::CellStatType type) const;
		void get_stat_by_real_cells(Stats::CellChrStatType stat, names_t &cell_barcodes, names_t &chromosome_names,
		                            counts_t &counts) const;

		const Cell &cell(size_t index) const;
		Cell &cell(size_t index);

		size_t has_exon_reads_num() const;
		size_t has_intron_reads_num() const;
		size_t has_not_annotated_reads_num() const;
		size_t real_cells_number() const;

		std::string merge_type() const;
		const StringIndexer& gene_indexer() const;
		const StringIndexer& umi_indexer() const;

		s_ul_hash_t umis_distribution() const;
	};
}
