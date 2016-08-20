#pragma once

#include <vector>

#include <Rcpp.h>

#include <Estimation/Stats.h>

namespace Estimation
{
	namespace Results
	{
		class CountMatrix
		{
		public:
			typedef std::vector<std::string> s_list_t;
			typedef std::vector<long> i_list_t;

		public:
			s_list_t cell_names;
			s_list_t gene_names;
			i_list_t counts;

			CountMatrix()
			{
			};

			CountMatrix(const s_list_t &cell_names, const s_list_t &gene_names, const i_list_t &counts)
					: cell_names(cell_names), gene_names(gene_names), counts(counts)
			{
			}
		};

		class IndropResult
		{
		public:
			typedef Stats::int_list_t int_list_t;

		public:
			CountMatrix cm;

			Stats::int_list_t ex_cells_chr_reads_counts;
			Stats::int_list_t nonex_cells_chr_reads_counts;
			Stats::str_list_t ex_cell_names;
			Stats::str_list_t nonex_cell_names;
			Stats::str_list_t chr_names;

			std::vector<double> reads_per_umi;
			Stats::int_list_t umig_covered;
			Stats::int_list_t merge_n;
			Stats::int_list_t reads_by_umig;
			Stats::str_list_t reads_by_umig_cbs;
			Stats::int_list_t exone_reads_by_cb;

			IndropResult()
			{ };

			IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
			             const int_list_t &umig_covered, bool not_filtered);

			virtual void save_rds(const std::string &filename) const;
			virtual Rcpp::List get_main_r_vec(const std::string &filename) const;
		};
	}
}