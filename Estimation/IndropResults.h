#pragma once

#include <vector>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

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
		ar & cm & non_exon_counts & non_exon_count_names & reads_per_umi & umig_covered & exon_counts &
		exon_count_names & merge_n;
	}

public:
	CountMatrix cm;
	std::vector<int> non_exon_counts;
	std::vector<std::string> non_exon_count_names;
	std::vector<int> exon_counts;
	std::vector<std::string> exon_count_names;
	std::vector<double> reads_per_umi;
	std::vector<int> umig_covered;
	std::vector<int> merge_n;

	IndropResult()
	{};

	IndropResult(const CountMatrix &cm, const std::vector<int> &non_exon_counts,
				 const std::vector<std::string> &non_exon_count_names, const std::vector<double> &reads_per_umi,
				 const std::vector<int> &umig_covered, const std::vector<int> &exon_counts,
				 const std::vector<std::string> &exon_count_names, const std::vector<int> &merge_n)
			: cm(cm)
			, non_exon_counts(non_exon_counts)
			, non_exon_count_names(non_exon_count_names)
			, reads_per_umi(reads_per_umi)
			, umig_covered(umig_covered)
			, exon_counts(exon_counts)
			, exon_count_names(exon_count_names)
			, merge_n(merge_n)
	{}
};
