#pragma once

#include "dropestr.h"
#include "umi_processing.h"

inline void hash_combine(std::size_t &seed) {}

template <typename T, typename... Rest>
inline void hash_combine(std::size_t &seed, const T &v, Rest... rest) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  hash_combine(seed, rest...);
}

struct pair_hash
{
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &p) const {
    size_t h;
    hash_combine(h, p.first, p.second);
    return h;
  }
};

class ClassifierData
{
private:
  using pairs_set = std::unordered_set<std::pair<std::string, std::string>, pair_hash>;

  class ClassifierRow
  {
  public:
    const std::string umi_small;
    const std::string umi_large;
    const int reads_per_umi_small;
    const int reads_per_umi_large;
    const char nucleotide_large;
    const int nucleotide_substitution;
    const int difference_position;
    const double umi_prob;
    const double log_quality_small;

    ClassifierRow(const std::string &umi_small, const std::string umi_large, int reads_per_umi_small,
                  int reads_per_umi_large, char nucleotide_large, int nucleotide_substitution,
                  int difference_position, double umi_prob, double log_quality_small);
  };

private:
  const sd_map_t _umi_probabilities_map;
  const sd_map_t _probability_normalizers;
  const bool _save_probs;

  std::vector<ClassifierRow> _data;
  pairs_set _presented_umi_pairs;

private:
  ClassifierRow fill_row(const UmiInfo &umi_small, const UmiInfo &umi_large, double umi_prob);

public:
  ClassifierData(const sd_map_t &umi_probabilities = sd_map_t(), const sd_map_t &probability_normalizers = sd_map_t());

  bool add_umis(const UmiInfo &umi1, const UmiInfo &umi2);
  Rcpp::DataFrame to_data_frame() const;
};
