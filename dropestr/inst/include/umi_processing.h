#pragma once

#include "Rcpp.h"
#include <vector>
#include <string>
#include <unordered_map>

class UmiInfo
{
public:
  typedef std::vector<double> quality_t;

private:
  std::string _umi;
  unsigned _reads_per_umi;
  quality_t _quality;

public:
  UmiInfo(const std::string &umi, const Rcpp::List &umi_data)
    : _umi(umi)
    , _reads_per_umi(Rcpp::as<unsigned>(umi_data[0]))
    , _quality(Rcpp::as<quality_t>(Rcpp::as<Rcpp::NumericVector>(umi_data[1]))) {}

  UmiInfo(const std::string &umi, const unsigned &reads_per_umi, const quality_t &quality)
    : _umi(umi)
    , _reads_per_umi(reads_per_umi)
    , _quality(quality) {}

  const std::string& umi() const {
    return this->_umi;
  }

  unsigned reads_per_umi() const {
    return this->_reads_per_umi;
  }

  const quality_t& quality() const {
    return this->_quality;
  }

  void merge(const UmiInfo &source) {
    if (this->_quality.size() != source.quality().size())
      Rcpp::stop("Quality vectors have different length");

    this->_reads_per_umi += source._reads_per_umi;

    for (int i = 0; i < this->_quality.size(); ++i) {
      this->_quality[i] = (this->_quality[i] * this->_reads_per_umi + source._quality[i] * source.reads_per_umi()) /
        (this->_reads_per_umi + source.reads_per_umi());
    }
  }

  Rcpp::List to_list() const {
    return Rcpp::List::create(this->_reads_per_umi, Rcpp::wrap(this->_quality));
  }
};

class UmisInfo
{
public:
  typedef std::unordered_map<std::string, UmiInfo> umis_info_t;

private:
  umis_info_t _info;

private:
  static umis_info_t parse_umis_info(const Rcpp::List &gene_info) {
    umis_info_t res;
    auto const umis = Rcpp::as<std::vector<std::string>>(Rcpp::as<Rcpp::StringVector>(gene_info.names()));

    for (int i = 0; i < gene_info.size(); ++i) {
      res.emplace(umis[i], UmiInfo(umis[i], Rcpp::as<Rcpp::List>(gene_info[i])));
    }

    return res;
  }

public:
  const umis_info_t& info() const {
    return this->_info;
  }

  UmisInfo()
    : _info()
  {}

  UmisInfo(const Rcpp::List &gene_info)
    : _info(parse_umis_info(gene_info))
  {}

  void add_umi(const UmiInfo &umi) {
    auto iter = this->_info.find(umi.umi());
    if (iter == this->_info.end()) {
      this->_info.emplace(umi.umi(), umi);
      return;
    }

    iter->second.merge(umi);
  }

  Rcpp::List to_list() const {
    Rcpp::List res(this->_info.size());
    Rcpp::StringVector umis(this->_info.size());

    size_t i = 0;
    for (auto umi_iter : this->_info) {
      res[i] = umi_iter.second.to_list();
      umis[i] = umi_iter.first;
      ++i;
    }

    res.attr("names") = umis;
    return res;
  }
};

class GeneInfo
{
public:
  typedef std::vector<unsigned> reads_per_umi_t;
  typedef std::vector<Rcpp::NumericVector> qualities_t;
  typedef std::vector<std::string> umis_t;

public:
  const umis_t umis;
  const reads_per_umi_t reads_per_umi;
  const qualities_t qualities;

private:
  static reads_per_umi_t parse_reads_per_umi(const Rcpp::List &gene_info) {
    reads_per_umi_t res(gene_info.size());
    for (int i = 0; i < gene_info.size(); ++i) {
      res[i] = Rcpp::as<Rcpp::List>(gene_info[i])[0];
    }

    return res;
  }

  static qualities_t parse_qualities(const Rcpp::List &gene_info) {
    qualities_t res(gene_info.size());
    for (int i = 0; i < gene_info.size(); ++i) {
      res[i] = Rcpp::as<Rcpp::List>(gene_info[i])[1];
    }

    return res;
  }

public:
  GeneInfo(const Rcpp::List &gene_info)
    : umis(Rcpp::as<umis_t>(Rcpp::as<Rcpp::StringVector>(gene_info.names())))
    , reads_per_umi(parse_reads_per_umi(gene_info))
    , qualities(parse_qualities(gene_info)) {}
};
