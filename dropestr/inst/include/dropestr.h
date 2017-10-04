#pragma once

#include <RcppEigen.h>

#include <progress.hpp>

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using s_set_t = std::unordered_set<std::string>;
using s_vec_t = std::vector<std::string>;
using ss_pair = std::pair<std::string, std::string>;
using si_map_t = std::unordered_map<std::string, int>;
using sd_map_t = std::unordered_map<std::string, double>;
using slst_map_t = std::unordered_map<std::string, Rcpp::List>;
using ssi_map_t = std::unordered_map<std::string, si_map_t>;
using umis_per_gene_t = std::unordered_map<std::string, ssi_map_t>;

const char NUCLEOTIDES[] = {'A', 'C', 'G', 'T'};
const int NUCLEOTIDES_NUM = 4;
const double EPS = 1e-8;

const si_map_t NUCL_PAIR_INDS = {
  std::make_pair("AC", 0),
  std::make_pair("AG", 1),
  std::make_pair("AT", 2),
  std::make_pair("CG", 3),
  std::make_pair("CT", 4),
  std::make_pair("GT", 5)
};

template<typename T>
s_vec_t as_s_vec(T vec) {
  return Rcpp::as<s_vec_t>(Rcpp::as<Rcpp::StringVector>(vec));
}

slst_map_t parseList(const Rcpp::List &lst);

si_map_t parseVector(const Rcpp::IntegerVector &vec);
sd_map_t parseVector(const Rcpp::NumericVector &vec);

Rcpp::NumericVector vpow(double base, const Rcpp::NumericVector& exp);
Rcpp::NumericVector vpow(const Rcpp::NumericVector& base, double exp);

si_map_t ValueCountsC(const s_vec_t &values);

Rcpp::List GetUmisDifference(const std::string &umi1, const std::string &umi2, int rpu1, int rpu2, bool force_neighbours = false, double umi_prob=-1);
