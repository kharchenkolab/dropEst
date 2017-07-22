#pragma once

#include <vector>

#include <RcppArmadillo.h>

#include <Estimation/Stats.h>

namespace Estimation
{
	class CellsDataContainer;

	class ResultsPrinter
	{
	private:
		typedef std::pair<std::string, unsigned long> si_pair_t;
		typedef std::vector<si_pair_t> s_counter_t;

		typedef std::vector<std::string> s_vec_t;
		typedef std::vector<int> l_vec_t;

	private:
		const bool text_output; //TODO: implement matrix output format
		const bool reads_output;

		static const size_t top_print_size = 10;

	private:
		static Rcpp::IntegerMatrix create_matrix(const s_vec_t &col_names, const s_vec_t &row_names,
		                                         const l_vec_t &counts);

		Rcpp::List get_saturation_analysis_info(const CellsDataContainer &container) const;
		Rcpp::DataFrame get_reads_per_chr_per_cell_info(Stats::CellStrStatType stat_type,
		                                                const CellsDataContainer &container,
		                                                const s_vec_t &cell_names) const;
		Rcpp::List get_reads_per_chr_per_cell_info(const CellsDataContainer &container, const s_vec_t &cell_names) const;
		s_vec_t get_filtered_cell_names(const CellsDataContainer &container) const;
		Rcpp::IntegerMatrix
		get_count_matrix(const CellsDataContainer &container, const s_counter_t &gene_counts) const;
		s_counter_t get_gene_counts_sorted(const CellsDataContainer &genes_container) const;
		Rcpp::NumericVector get_mean_reads_per_umi(const CellsDataContainer &container) const;
		Rcpp::List get_reads_per_umi_per_cell(const CellsDataContainer &container) const;
		Rcpp::List get_merge_targets(const CellsDataContainer &container) const;

	public:
		ResultsPrinter(bool text_output, bool reads_output);

		void save_results(const CellsDataContainer &container, const std::string &filename) const;
	};
}