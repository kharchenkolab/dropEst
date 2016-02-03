#pragma once

#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

#ifdef R_LIBS
#include <Rcpp.h>
#endif

class Stats;

class CountMatrix
{
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int /* file_version */)
	{
		ar & cell_names & gene_names & counts;
	}

public:
	std::vector<std::string> cell_names;
	std::vector<std::string> gene_names;
	std::vector<int> counts;

	CountMatrix()
	{};

	CountMatrix(const std::vector<std::string> &cell_names, const std::vector<std::string> &gene_names,
				const std::vector<int> &counts)
			: cell_names(cell_names)
			, gene_names(gene_names)
			, counts(counts)
	{}
};

class IndropResult
{
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int /* file_version */)
	{
		ar & this->cm & this->non_exon_chr_counts & this->non_exon_chr_count_names & this->reads_per_umi &
				this->umig_covered & this->exon_chr_counts & this->exon_chr_count_names & this->merge_n;
	}

public:
	CountMatrix cm;
	std::vector<int> non_exon_chr_counts;
	std::vector<std::string> non_exon_chr_count_names;
	std::vector<int> exon_chr_counts;
	std::vector<std::string> exon_chr_count_names;

	std::vector<int> non_exon_cell_counts;
	std::vector<std::string> non_exon_cell_count_tags;
	std::vector<int> exon_cell_counts;
	std::vector<std::string> exon_cell_count_tags;
	
	std::vector<double> reads_per_umi;
	std::vector<int> umig_covered;
	std::vector<int> merge_n;

	IndropResult() = default;

	IndropResult(const CountMatrix &cm, const std::vector<int> &non_exon_chr_counts,
				 const std::vector<std::string> &non_exon_chr_count_names, const std::vector<int> &exon_chr_counts,
				 const std::vector<std::string> &exon_chr_count_names, const std::vector<double> &reads_per_umi,
				 const std::vector<int> &umig_covered, const std::vector<int> &merge_n);

	IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
				 const std::vector<int> &umig_covered);

#ifdef R_LIBS
	Rcpp::List get_r_table(const std::string &fname) const;
#endif
};
