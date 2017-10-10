#pragma once

#include <Rcpp.h>
#include <progress.hpp>

#include <omp.h>
// #include <gperftools/profiler.h>

#include <set>
#include <unordered_map>
#include <unordered_set>
using namespace Rcpp;

using s_set_t = std::unordered_set<std::string>;
using s_vec_t = std::vector<std::string>;
using ss_pair = std::pair<std::string, std::string>;
using si_map_t = std::unordered_map<std::string, int>;
using sd_map_t = std::unordered_map<std::string, double>;
using ssi_map_t = std::unordered_map<std::string, si_map_t>;
using umis_per_gene_t = std::unordered_map<std::string, ssi_map_t>;
using neighbours_map_t = std::unordered_map<std::string, s_vec_t>;

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
  return as<s_vec_t>(as<StringVector>(vec));
}

si_map_t parseVector(const IntegerVector &vec) {
  si_map_t res;
  StringVector names = as<StringVector>(vec.names());
  for (int i = 0; i < vec.size(); ++i) {
    res.emplace(std::string(names[i]), vec[i]);
  }

  return res;
}

sd_map_t parseVector(const NumericVector &vec) {
  sd_map_t res;
  StringVector names = as<StringVector>(vec.names());
  for (int i = 0; i < vec.size(); ++i) {
    res.emplace(std::string(names[i]), vec[i]);
  }

  return res;
}

NumericVector vpow(double base, const NumericVector& exp) {
  NumericVector res(exp.size());
  for (int i = 0; i < res.size(); ++i) {
    res[i] = std::pow(base, exp[i]);
  }
  return res;
}

NumericVector vpow(const NumericVector& base, double exp) {
  NumericVector res(base.size());
  for (int i = 0; i < res.size(); ++i) {
    res[i] = std::pow(base[i], exp);
  }
  return res;
}

// [[Rcpp::export]]
si_map_t ValueCountsC(const s_vec_t &values) {
  si_map_t res;
  for (auto const value : values) {
    res[value]++;
  }

  return res;
}

// [[Rcpp::export]]
List GetUmisDifference(const std::string &umi1, const std::string &umi2, int rpu1, int rpu2, bool force_neighbours = false, double umi_prob=-1) {
  if (umi1.length() != umi2.length())
    stop("UMIs must have the same length");

  int diff_pos = -1, diff_num = 0;
  std::string nuc_diff = "NN";
  for (int i = 0; i < umi1.length(); ++i) {
    char n1 = umi1[i], n2 = umi2[i];
    if (n1 == n2)
      continue;

    diff_num++;
    if (diff_num != 1)
      continue;

    diff_pos = i;
    nuc_diff[0] = std::min(n1, n2);
    nuc_diff[1] = std::max(n1, n2);

    if (force_neighbours)
      break;
  }

  if (diff_num > 1) {
    diff_pos = runif(1, 0, umi1.length())[0];
    char n1, n2;
    do {
      n1 = NUCLEOTIDES[int(runif(1, 0, NUCLEOTIDES_NUM)[0])];
      n2 = NUCLEOTIDES[int(runif(1, 0, NUCLEOTIDES_NUM)[0])];
    }
    while (n1 == n2);

    nuc_diff[0] = std::min(n1, n2);
    nuc_diff[1] = std::max(n1, n2);
  }

  List res = List::create(_["Position"]=diff_pos, _["Nucleotides"]=NUCL_PAIR_INDS.at(nuc_diff), _["ED"]=diff_num,
                       _["MinRpU"]=std::min(rpu1, rpu2), _["MaxRpU"]=std::max(rpu1, rpu2));

  if (umi_prob > 0) {
    res["UmiProb"] = umi_prob;
  }

  return res;
}