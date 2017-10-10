#pragma once

#include <random>
#include "CommonFuncs.hpp"

int GetGeneEstimatedValueClassic(int value, int umis_number);

// unsigned hammingDistance(const std::string &s1, const std::string &s2, bool only_first=false) {
//   if (s1.size() != s2.size())
//     throw std::runtime_error("Strings must have equal length");
//
//   int dist = 0;
//   for (int i = 0; i < s1.size(); ++i) {
//     if (s1[i] != s2[i]) {
//       dist++;
//       if (only_first && dist > 1)
//         break;
//     }
//   }
//
//   return dist;
// }
//
// // [[Rcpp::export]]
// double GetBootstrapAdjacentUmisProbability(const std::vector<double> &umi_probabilities, const std::vector<std::string> &umis,
//                                            unsigned long repeats_num = 10000, int seed = -1) {
//   if (seed == -1) {
//     seed = time(nullptr);
//   }
//
//   srand(seed);
//
//   std::discrete_distribution<int> dist(umi_probabilities.begin(), umi_probabilities.end());
//   std::mt19937 gen;
//
//   unsigned long adjacent_num = 0;
//   for (unsigned long i = 0; i < repeats_num; ++i) {
//     int i1 = dist(gen), i2 = dist(gen);
//     unsigned ed = hammingDistance(umis[i1], umis[i2], true);
//     if (ed == 1) {
//       adjacent_num++;
//     }
//   }
//
//   return double(adjacent_num) / repeats_num;
// }

// [[Rcpp::export]]
List GetTrimCollisionsNum(const List &rpu_per_cell, int trim_length) {
  List res(rpu_per_cell.size());

  auto const &cell_names = as<s_vec_t>(as<StringVector>(rpu_per_cell.names()));
  for (int cell_ind = 0; cell_ind < rpu_per_cell.size(); ++cell_ind) {
    const auto &rpu_per_gene = as<List>(rpu_per_cell[cell_ind]);
    IntegerVector collisions_per_gene(rpu_per_gene.size());

    for (int gene_ind = 0; gene_ind < rpu_per_gene.size(); ++gene_ind) {
      auto const &reads_per_umi = as<IntegerVector>(rpu_per_gene[gene_ind]);
      auto const &umis = as<s_vec_t>(as<StringVector>(reads_per_umi.names()));
      s_set_t trimmed_rpus;

      int collisions_num = 0;
      for (int umi_ind = 0; umi_ind < reads_per_umi.size(); ++umi_ind) {
        auto ins = trimmed_rpus.insert(umis[umi_ind].substr(0, trim_length));
        if (!ins.second) {
          collisions_num++;
        }
      }

      collisions_per_gene[gene_ind] = collisions_num;
    }

    collisions_per_gene.attr("names") = as<StringVector>(rpu_per_gene.names());
    res[cell_ind] = collisions_per_gene;
  }

  res.attr("names") = cell_names;
  return res;
}

// [[Rcpp::export]]
double GetBootstrapUmisMeanNum(const std::vector<double> &umi_probabilities, size_t size, unsigned repeats_num, int seed = -1) {
  if (seed == -1) {
    seed = time(nullptr);
  }

  srand(seed);

  std::discrete_distribution<int> dist(umi_probabilities.begin(), umi_probabilities.end());
  std::mt19937 gen;

  unsigned long sum = 0;
  for (unsigned i = 0; i < repeats_num; ++i) {
    std::vector<int> uniq_vals(size);
    for (int try_num = 0; try_num < size; ++try_num) {
      uniq_vals[try_num] = dist(gen);
    }
    std::sort(uniq_vals.begin(), uniq_vals.end());
    int diff_num = std::distance(uniq_vals.begin(), std::unique(uniq_vals.begin(), uniq_vals.end()));
    sum += diff_num;
  }

  return double(sum) / repeats_num;
}

int roundInt(int val, int order) {
  return std::round(val / double(order)) * order;
}

// [[Rcpp::export]]
int GetGeneEstimatedValueClassic(int value, int umis_number) {
  if (value == umis_number) {
    return 2 * GetGeneEstimatedValueClassic(value - 1, umis_number) - GetGeneEstimatedValueClassic(value - 2, umis_number);
  }
  return int(std::round(-std::log(1 - value / double(umis_number)) * umis_number));
}

// [[Rcpp::export]]
int GetGeneEstimatedValue(int value, const std::vector<int> &real_sizes, const std::vector<double> &estimated_sizes) {
  if (estimated_sizes.empty() || real_sizes.empty() || real_sizes.size() != estimated_sizes.size())
    stop("Bad input arrays: " + std::to_string(real_sizes.size()) + ", " + std::to_string(estimated_sizes.size()));

  if (estimated_sizes.back() == value)
    return real_sizes.back();

  auto est_iter = std::upper_bound(estimated_sizes.begin(), estimated_sizes.end(), double(value));
  if (est_iter == estimated_sizes.end())
    stop("Too large value of gene expression: " + std::to_string(value));

  int index = std::distance(estimated_sizes.begin(), est_iter);
  if (index == 0)
    stop("Undexpected expression value: " + std::to_string(value));

  int lower_real = real_sizes[index - 1], upper_real = real_sizes[index];
  double lower_est = estimated_sizes[index - 1], upper_est = estimated_sizes[index];

  return lower_real + (upper_real - lower_real) * (value - lower_est) / (upper_est - lower_est); // Linear interpolation
}