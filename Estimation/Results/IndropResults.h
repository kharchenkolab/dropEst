#pragma once

#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>
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

		private:
			friend class boost::serialization::access;

			template<class Archive>
			void serialize(Archive &ar, const unsigned int /* file_version */)
			{
				ar & cell_names & gene_names & counts;
			}

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
			friend class boost::serialization::access;

		protected:
			template<class Archive>
			void serialize(Archive &ar, const unsigned int /* file_version */)
			{
				ar & this->cm & this->ex_cells_chr_reads_counts & this->nonex_cells_chr_reads_counts &
				this->ex_cell_names &
				this->nonex_cell_names & this->chr_names & this->reads_per_umi & this->umig_covered & this->merge_n &
				this->reads_by_umig;
			}

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
			int_list_t umig_covered;
			int_list_t merge_n;
			int_list_t reads_by_umig;
			int_list_t exone_reads_by_cb;

			IndropResult()
			{ };

			IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
			             const int_list_t &umig_covered, bool not_filtered);

#ifdef R_LIBS
			virtual void save_rds(const std::string &filename) const;
			virtual Rcpp::List get_main_r_vec(const std::string &filename) const;
#endif
		};
	}
}