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
		typedef std::vector<std::string> s_vec_t;
		typedef std::vector<int> i_vec_t;
		typedef std::vector<arma::uword> arma_uvec;

	private:
		const bool write_matrix;
		const bool reads_output;

		static const size_t top_print_size = 10;

	private:
		static Rcpp::IntegerMatrix create_matrix(const s_vec_t &col_names, const s_vec_t &row_names,
		                                         const i_vec_t &counts);
		static arma::sp_umat create_matrix(const arma_uvec &row_numbers, const arma_uvec &column_numbers,
		                                   const arma_uvec &counts, size_t total_rows, size_t total_cols);

		Rcpp::List get_saturation_analysis_info(const CellsDataContainer &container) const;
		Rcpp::DataFrame get_reads_per_chr_per_cell_info(Stats::CellStrStatType stat_type,
		                                                const CellsDataContainer &container,
		                                                const s_vec_t &cell_names) const;
		Rcpp::List get_reads_per_chr_per_cell_info(const CellsDataContainer &container, const s_vec_t &cell_names) const;
		s_vec_t get_filtered_cell_names(const CellsDataContainer &container) const;
		s_vec_t get_real_cell_names(const CellsDataContainer &container) const;
		Rcpp::List get_count_matrix(const CellsDataContainer &container, bool filtered) const;
		void trace_gene_counts(const CellsDataContainer &genes_container) const;
		Rcpp::NumericVector get_mean_reads_per_umi(const CellsDataContainer &container) const;
		Rcpp::List get_reads_per_umi_per_cell(const CellsDataContainer &container) const;
		Rcpp::List get_merge_targets(const CellsDataContainer &container) const;

		arma::sp_umat get_count_matrix_filtered(const CellsDataContainer &container, s_vec_t &gene_names,
		                                        s_vec_t &cell_names) const;

		arma::sp_umat get_count_matrix_raw(const CellsDataContainer &container, s_vec_t &gene_names,
		                                        s_vec_t &cell_names) const;

	public:
		ResultsPrinter(bool write_matrix, bool reads_output);

		void save_results(const CellsDataContainer &container, const std::string &filename) const;

		Rcpp::IntegerVector get_requested_umis_per_cb(const CellsDataContainer &container) const;
	};
}
