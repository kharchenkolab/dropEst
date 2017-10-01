#pragma once

#include "Rcpp.h"

class UmiInfo
{
public:
	typedef std::vector<int> quality_t;
public:
	const std::string umi;
	const unsigned reads_per_umi;
	const quality_t quality;

	UmiInfo(const std::string &umi, const Rcpp::List &umi_data)
			: umi(umi)
			, reads_per_umi(Rcpp::as<unsigned>(umi_data[0]))
			, quality(Rcpp::as<quality_t>(Rcpp::as<Rcpp::IntegerVector>(umi_data[1])))
	{}

	UmiInfo(const std::string &umi, const unsigned &reads_per_umi, const std::vector<int> &quality)
		: umi(umi)
		, reads_per_umi(reads_per_umi)
		, quality(quality)
	{}
};

class UmisInfo
{
public:
	typedef std::vector<UmiInfo> umis_info_vec_t;
private:
	static umis_info_vec_t parse_umis_info(const Rcpp::List &gene_info)
	{
		umis_info_vec_t res;
		auto const umis = Rcpp::as<std::vector<std::string>>(Rcpp::as<Rcpp::StringVector>(gene_info.names()));

		for (int i = 0; i < gene_info.size(); ++i)
		{
			res.push_back(UmiInfo(umis[i], Rcpp::as<Rcpp::List>(gene_info[i])));
		}

		return res;
	}

public:
	const umis_info_vec_t info;
	UmisInfo(const Rcpp::List &gene_info)
		: info(parse_umis_info(gene_info))
	{}
};

class GeneInfo
{
public:
	typedef std::vector<unsigned> reads_per_umi_t;
	typedef std::vector<Rcpp::IntegerVector> qualities_t;

public:
	const Rcpp::StringVector umis;
	const reads_per_umi_t reads_per_umi;
	const qualities_t qualities;

private:
	static reads_per_umi_t parse_reads_per_umi(const Rcpp::List &gene_info)
	{
		reads_per_umi_t res(gene_info.size());
		for (int i = 0; i < gene_info.size(); ++i)
		{
			res[i] = Rcpp::as<Rcpp::List>(gene_info[i])[0];
		}

		return res;
	}

	static qualities_t parse_qualities(const Rcpp::List &gene_info)
	{
		qualities_t res(gene_info.size());
		for (int i = 0; i < gene_info.size(); ++i)
		{
			res[i] = Rcpp::as<Rcpp::List>(gene_info[i])[1];
		}

		return res;
	}

public:
	GeneInfo(const Rcpp::List& gene_info)
		: umis(Rcpp::as<Rcpp::StringVector>(gene_info.names()))
		, reads_per_umi(parse_reads_per_umi(gene_info))
		, qualities(parse_qualities(gene_info))
	{}
};