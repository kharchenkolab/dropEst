#pragma once

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

#include "CellsDataContainer.h"

namespace Estimation
{
	namespace Results
	{
		class CountMatrix;
		class IndropResult;
		class BadCellsStats;
	}

	namespace Merge
	{
		namespace UMIs
		{
			class MergeUMIsStrategySimple;
		}
	}

	class Estimator
	{
	private:
		typedef boost::unordered_map<std::string, int> SIHM;
		typedef boost::unordered_set<std::string> s_set;

		typedef std::vector<std::pair<std::string, int> > s_counter_t;
		typedef std::vector<double> doubles_t;
		typedef CellsDataContainer::names_t names_t;
		typedef std::vector<int> i_list_t;
		typedef std::vector<long> l_list_t;
		typedef std::shared_ptr<Merge::MergeStrategyAbstract> merge_strat_ptr;
		typedef std::shared_ptr<Merge::UMIs::MergeUMIsStrategySimple> umi_merge_strat_ptr;

	public:
		typedef boost::unordered_map<std::string, boost::unordered_map<std::string, unsigned>> ss_u_hash_t; //Can't use long because of RInside (see BadCellsStats)

	private:
		static const size_t top_print_size = 10;

		const merge_strat_ptr merge_strategy;
		const umi_merge_strat_ptr umi_merge_strategy;

	public:
		Estimator(const merge_strat_ptr &merge_strategy, const umi_merge_strat_ptr &umi_merge_strategy);

		Results::IndropResult get_results(const CellsDataContainer &container, bool not_filtered, bool reads_output);
		Results::BadCellsStats get_bad_cells_results(const CellsDataContainer &container);

		CellsDataContainer get_cells_container(const names_t &files, bool bam_output, bool filled_bam, int max_cells_num,
	                                           const std::string &reads_params_names_str, const std::string &gtf_filename,
				                               const std::vector<CellsDataContainer::Mark> &gene_match_levels);

	private:
		names_t get_filtered_cell_names(const CellsDataContainer &genes_container) const;

		i_list_t get_count_matrix(const names_t &gene_names, const CellsDataContainer &genes_container, bool reads_output) const;

		names_t get_gene_names_sorted(const CellsDataContainer &genes_container) const;

		doubles_t get_reads_per_umis(const CellsDataContainer &genes_container) const;

		l_list_t get_umig_coverage(const CellsDataContainer &genes_container) const;

		Results::IndropResult get_indrop_results(const Results::CountMatrix &cm, const CellsDataContainer &genes_conteiner,
												 bool not_filtered) const;

		std::string get_genes_top_verbose(const s_counter_t &genes_count) const;
	};
}
