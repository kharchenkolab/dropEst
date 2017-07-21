#include "dropestr.h"

using namespace Rcpp;

inline void hash_combine(std::size_t& seed) { }

template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  hash_combine(seed, rest...);
}

struct pair_hash {
  template <class T1, class T2>
  std::size_t operator () (const std::pair<T1,T2> &p) const {
    size_t h;
    hash_combine(h, p.first, p.second);
    return h;
  }
};

// [[Rcpp::export]]
NumericVector PredictNBC(const List &clf, const List &predict_data) {
  std::vector<NumericVector> probs(2);
  for (int i = 0; i < 2; ++i) {
    auto const &cur_clf = as<List>(clf[i]);

    NumericVector nucl_prob = as<NumericVector>(cur_clf["Nucleotides"])[as<IntegerVector>(predict_data["Nucleotides"])];
    NumericVector position_prob = as<NumericVector>(cur_clf["Position"])[as<IntegerVector>(predict_data["Position"])];

    NumericVector min_rpu_prob = vpow(as<double>(cur_clf["MinRpU"]), as<NumericVector>(predict_data["MinRpU"]));
    NumericVector max_rpu_prob = vpow(as<double>(cur_clf["MaxRpU"]), as<NumericVector>(predict_data["MaxRpU"]));

    probs[i] = nucl_prob * position_prob * min_rpu_prob * max_rpu_prob;
  }

  return probs[1] / (probs[1] + probs[0]);
}

// [[Rcpp::export]]
NumericVector PredictLeftPart(const List &neg_clf, const std::vector<double> &rpu_prob, const DataFrame &predict_data, int gene_size) {
  NumericVector nucl_prob = log(as<NumericVector>(as<NumericVector>(neg_clf["Nucleotides"])[as<IntegerVector>(predict_data["Nucleotides"])]));
  NumericVector position_prob = log(as<NumericVector>(as<NumericVector>(neg_clf["Position"])[as<IntegerVector>(predict_data["Position"])]));

  NumericVector min_rpu_prob_neg = std::log(as<double>(neg_clf["MinRpU"])) * as<NumericVector>(predict_data["MinRpU"]);
  NumericVector max_rpu_prob_neg = std::log(as<double>(neg_clf["MaxRpU"])) * as<NumericVector>(predict_data["MaxRpU"]);

  NumericVector min_rpu_prob_pos(max_rpu_prob_neg.size());
  NumericVector max_rpu_prob_pos(max_rpu_prob_neg.size());

  for (int i = 0; i < max_rpu_prob_neg.size(); ++i) {
    min_rpu_prob_pos[i] = std::log(rpu_prob.at(as<IntegerVector>(predict_data["MinRpU"])[i] - 1));
    max_rpu_prob_pos[i] = std::log(rpu_prob.at(as<IntegerVector>(predict_data["MaxRpU"])[i] - 1));
  }

  NumericVector umi_prob_pos = log(NumericVector(1 - vpow(1 - as<NumericVector>(predict_data["UmiProb"]), gene_size)));

  return exp((nucl_prob + position_prob + min_rpu_prob_neg + max_rpu_prob_neg) - (min_rpu_prob_pos + max_rpu_prob_pos + umi_prob_pos));
}

// [[Rcpp::export]]
std::vector<int> ArrangePredictions(const std::vector<int> &target_umi_factors, const std::vector<double> &probs) {
  if (target_umi_factors.size() != probs.size())
    stop("Size must be equal");

  std::vector<int> res(target_umi_factors.size());
  std::iota(res.begin(), res.end(), 1);

  std::sort(res.begin(), res.end(), [target_umi_factors, probs](int i, int j){
    int f1 = target_umi_factors[i - 1], f2 = target_umi_factors[j - 1];
    if (f1 != f2)
      return f1 < f2;

    return probs[i - 1] < probs[j - 1];
  });

  return res;
}

// [[Rcpp::export]]
std::vector<bool> FilterPredictions(const s_vec_t &not_filtered_umis, const s_vec_t &base_umis, const s_vec_t &target_umis) {
  s_set_t not_filtered_umis_set(not_filtered_umis.begin(), not_filtered_umis.end());
  std::vector<bool> res(base_umis.size(), false);
  for (int i = 0; i < res.size(); ++i) {
    if (not_filtered_umis_set.find(base_umis[i]) == not_filtered_umis_set.end())
      continue;

    if (not_filtered_umis_set.find(target_umis[i]) == not_filtered_umis_set.end())
      continue;

    res[i] = true;
  }

  return res;
}

// [[Rcpp::export]]
List PrepareClassifierData(const IntegerVector &reads_per_umi, const List &neighborhood, const NumericVector &umi_probabilities = NumericVector()) {
  auto rpu_map = parseVector(reads_per_umi);
  sd_map_t umi_probabs_map;
  if (umi_probabilities.size() > 0) {
    umi_probabs_map = parseVector(umi_probabilities);
  }

  int res_size = 0;
  for (auto const &neighbours : neighborhood) {
    res_size += as<StringVector>(neighbours).size();
  }

  using pairs_set = std::unordered_set<std::pair<std::string, std::string>, pair_hash>;
  pairs_set umi_pairs_set;

  List differences(res_size);
  std::vector<pairs_set::const_iterator> umi_pairs;

  auto const &umis = as<s_vec_t>(neighborhood.names());
  for (int i = 0; i < neighborhood.size(); ++i) {
    const std::string &umi1 = umis[i];
    int rpu1 = rpu_map.at(umi1);
    for (auto const &umi2 : as<s_vec_t>(as<StringVector>(neighborhood[i]))) {
      int rpu2 = rpu_map.at(umi2);

      ss_pair cur_umi_pair;
      if (rpu1 <= rpu2) {
        cur_umi_pair = std::make_pair(umi1, umi2);
      }
      else {
        cur_umi_pair = std::make_pair(umi2, umi1);
      }

      auto iter = umi_pairs_set.insert(cur_umi_pair);
      if (!iter.second)
        continue;

      double umi_prob = -1;
      if (!umi_probabs_map.empty()) {
        umi_prob = umi_probabs_map.at(cur_umi_pair.first);
      }

      differences[umi_pairs.size()] = GetUmisDifference(cur_umi_pair.first, cur_umi_pair.second, rpu1, rpu2, true, umi_prob);
      umi_pairs.push_back(iter.first);
    }
  }

  StringVector umi_matrix(Dimension(umi_pairs.size(), 2));

  for (int i = 0; i < umi_pairs.size(); ++i) {
    umi_matrix(i, 0) = umi_pairs[i]->first;
    umi_matrix(i, 1) = umi_pairs[i]->second;
  }

  colnames(umi_matrix) = StringVector::create("Base", "Target");

  if (umi_pairs.size() == 0) {
    differences = List();
  }
  else {
    differences = differences[Range(0, umi_pairs.size() - 1)];
  }

  return List::create(_["umi.pairs"]=DataFrame(umi_matrix), _["differences"]=differences);
}
