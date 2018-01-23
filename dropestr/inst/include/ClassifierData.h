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
    const double log_quality_small;

    ClassifierRow(const std::string &umi_small, const std::string umi_large, int reads_per_umi_small,
                  int reads_per_umi_large, double log_quality_small);
  };

private:
  std::vector<ClassifierRow> _data;
  pairs_set _presented_umi_pairs;

private:
  ClassifierRow fill_row(const UmiInfo &umi_small, const UmiInfo &umi_large);

public:
  bool add_umis(const UmiInfo &umi1, const UmiInfo &umi2);
  Rcpp::DataFrame to_data_frame() const;
};
