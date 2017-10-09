#include <umi_processing.h>
#include "dropestr.h"

using namespace Rcpp;

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
  class ClassifierRow
  {
  public:
    const std::string umi_small;
    const std::string umi_large;
    const int reads_per_umi_small;
    const int reads_per_umi_large;
    const int edit_distance;
    const int nucleotide_substitution;
    const int difference_position;
    const double umi_prob;
    const double log_quality_small;

    ClassifierRow(const std::string &umi_small, const std::string umi_large, int reads_per_umi_small,
                  int reads_per_umi_large, int edit_distance, int nucleotide_substitution, int difference_position,
                  double umi_prob, double log_quality_small)
      : umi_small(umi_small)
      , umi_large(umi_large)
      , reads_per_umi_small(reads_per_umi_small)
      , reads_per_umi_large(reads_per_umi_large)
      , edit_distance(edit_distance)
      , nucleotide_substitution(nucleotide_substitution)
      , difference_position(difference_position)
      , umi_prob(umi_prob)
      , log_quality_small(log_quality_small) {}
  };

  using pairs_set = std::unordered_set<std::pair<std::string, std::string>, pair_hash>;

private:
  const sd_map_t _umi_probabilities_map;
  const bool _save_probs;
  const bool _force_neighbours;
  std::vector<ClassifierRow> _data;
  pairs_set _presented_umi_pairs;

private:
  ClassifierRow fill_row(const UmiInfo &umi_small, const UmiInfo &umi_large, double umi_prob) {
    if (umi_small.umi().length() != umi_large.umi().length())
      stop("UMIs must have the same length");

    int diff_pos = -1, diff_num = 0;
    std::string nuc_diff = "NN";
    for (int i = 0; i < umi_small.umi().length(); ++i) {
      char n1 = umi_small.umi()[i], n2 = umi_large.umi()[i];
      if (n1 == n2)
        continue;

      diff_num++;
      if (diff_num != 1)
        continue;

      diff_pos = i;
      nuc_diff[0] = std::min(n1, n2);
      nuc_diff[1] = std::max(n1, n2);

      if (this->_force_neighbours)
        break;
    }

    if (diff_num > 1) {
      diff_pos = runif(1, 0, umi_small.umi().length())[0];
      char n1, n2;
      do {
        n1 = NUCLEOTIDES[int(runif(1, 0, NUCLEOTIDES_NUM)[0])];
        n2 = NUCLEOTIDES[int(runif(1, 0, NUCLEOTIDES_NUM)[0])];
      } while (n1 == n2);

      nuc_diff[0] = std::min(n1, n2);
      nuc_diff[1] = std::max(n1, n2);
    }

    if (diff_pos == -1) {
      warning("Equal UMIs were compared");
    }

    if (diff_num == 1) {
      double quality_small = diff_pos == -1 ? -1 : umi_small.quality()[diff_pos];
      return ClassifierRow(umi_small.umi(), umi_large.umi(), umi_small.reads_per_umi(), umi_large.reads_per_umi(),
                           diff_num, NUCL_PAIR_INDS.at(nuc_diff), diff_pos, umi_prob, quality_small);
    }

    return ClassifierRow(umi_small.umi(), umi_large.umi(), umi_small.reads_per_umi(), umi_large.reads_per_umi(),
                         diff_num, NA_INTEGER, NA_INTEGER, umi_prob, NA_REAL);
  }

public:
  bool add_umis(const UmiInfo &umi1, const UmiInfo &umi2) {
    const UmiInfo &umi_small = umi1.reads_per_umi() <= umi2.reads_per_umi() ? umi1 : umi2;
    const UmiInfo &umi_large = umi1.reads_per_umi() <= umi2.reads_per_umi() ? umi2 : umi1;

    if (!this->_presented_umi_pairs.insert(std::make_pair(umi_small.umi(), umi_large.umi())).second)
      return false;

    double umi_prob = -1;
    if (!this->_umi_probabilities_map.empty()) {
      umi_prob = this->_umi_probabilities_map.at(umi_small.umi());
    }
    this->_data.push_back(this->fill_row(umi_small, umi_large, umi_prob));
    return true;
  }

  ClassifierData(bool force_neighbours, const sd_map_t &umi_probabilities_map = sd_map_t())
    : _umi_probabilities_map(umi_probabilities_map)
    , _save_probs(!_umi_probabilities_map.empty())
    , _force_neighbours(force_neighbours)
  {}

  DataFrame to_data_frame() const {
    if (this->_data.empty())
      return DataFrame();

    StringVector umi_small(this->_data.size());
    StringVector umi_large(this->_data.size());
    IntegerVector reads_per_umi_small(this->_data.size());
    IntegerVector reads_per_umi_large(this->_data.size());
    IntegerVector edit_distance(this->_data.size());
    IntegerVector nucleotide_substitution(this->_data.size());
    IntegerVector difference_position(this->_data.size());
    NumericVector umi_prob(this->_data.size());
    NumericVector mean_log_quality_small(this->_data.size());

    for (size_t i = 0; i < this->_data.size(); ++i) {
      auto const &umi_info = this->_data[i];
      umi_small[i] = umi_info.umi_small;
      umi_large[i] = umi_info.umi_large;
      reads_per_umi_small[i] = umi_info.reads_per_umi_small;
      reads_per_umi_large[i] = umi_info.reads_per_umi_large;
      edit_distance[i] = umi_info.edit_distance;
      nucleotide_substitution[i] = umi_info.nucleotide_substitution;
      difference_position[i] = umi_info.difference_position;
      umi_prob[i] = umi_info.umi_prob;
      mean_log_quality_small[i] = umi_info.log_quality_small;
    }

    auto res = DataFrame::create(_["Base"] = umi_small, _["Target"] = umi_large, _["Position"] = difference_position,
                                 _["Nucleotides"] = nucleotide_substitution, _["ED"] = edit_distance,
                                 _["MinRpU"] = reads_per_umi_small, _["MaxRpU"] = reads_per_umi_large,
                                 _["Quality"] = mean_log_quality_small);

    if (this->_save_probs) {
      res["UmiProb"] = umi_prob;
    }

    return res;
  }
};

//' @export
// [[Rcpp::export]]
IntegerVector Quantize(const NumericVector &values, const NumericVector &quant_borders) {
  const double EPS = 1e-7;
  IntegerVector quants(values.size());

  for (int val_id = 0; val_id < values.size(); ++val_id) {
    const double value = values[val_id];
    for (int border_id = 0; border_id < quant_borders.size(); ++border_id) {
      if (value < quant_borders[border_id] + EPS || border_id == quant_borders.size() - 1) {
        quants[val_id] = border_id;
        break;
      }
    }
  }

  return quants;
}

// [[Rcpp::export]]
NumericVector PredictLeftPart(const List &classifier, const DataFrame &predict_data, int gene_size) {
  using NV = NumericVector;
  using IV = IntegerVector;

  List error_distributions = as<List>(classifier["Negative"]);
  List common_distributions = as<List>(classifier["Common"]);

  NV nucl_prob_err = log(as<NV>(as<NV>(error_distributions["Nucleotides"])[as<IV>(predict_data["Nucleotides"])]));
  NV position_prob_err = log(as<NV>(as<NV>(error_distributions["Position"])[as<IV>(predict_data["Position"])]));

  NV min_rpu_prob_err = std::log(as<double>(error_distributions["MinRpU"])) * as<NV>(predict_data["MinRpU"]);
  NV max_rpu_prob_err = std::log(as<double>(error_distributions["MaxRpU"])) * as<NV>(predict_data["MaxRpU"]);

  IV quantized_quality = Quantize(as<NV>(predict_data["Quality"]), as<NV>(classifier["QualityQuantBorders"]));

  NV rpu_prob = as<NV>(common_distributions["RpuProbs"]);

  NV min_rpu_prob(max_rpu_prob_err.size());
  NV max_rpu_prob(max_rpu_prob_err.size());
  NV quality_prob(max_rpu_prob_err.size());
  NV quality_prob_err(max_rpu_prob_err.size());

  for (int i = 0; i < max_rpu_prob_err.size(); ++i) {
    min_rpu_prob[i] = std::log(rpu_prob.at(as<IV>(predict_data["MinRpU"])[i] - 1));
    max_rpu_prob[i] = std::log(rpu_prob.at(as<IV>(predict_data["MaxRpU"])[i] - 1));

    int current_quality_quant = quantized_quality[i];
    quality_prob[i] = std::log(as<NV>(common_distributions["Quality"]).at(current_quality_quant));
    quality_prob_err[i] = std::log(as<NV>(error_distributions["Quality"]).at(current_quality_quant));
  }

  NV umi_prob_pos = log(NV(1 - vpow(1 - as<NV>(predict_data["UmiProb"]), gene_size)));

  return exp((nucl_prob_err + position_prob_err + min_rpu_prob_err + max_rpu_prob_err) - //  + quality_prob_err
             (umi_prob_pos + min_rpu_prob + max_rpu_prob)); //  + quality_prob
}

// [[Rcpp::export]]
std::vector<int> ArrangePredictions(const std::vector<int> &target_umi_factors, const std::vector<double> &probs) {
  if (target_umi_factors.size() != probs.size())
    stop("Size must be equal");

  std::vector<int> res(target_umi_factors.size());
  std::iota(res.begin(), res.end(), 1);

  std::sort(res.begin(), res.end(), [target_umi_factors, probs](int i, int j) {
    int f1 = target_umi_factors[i - 1], f2 = target_umi_factors[j - 1];
    if (f1 != f2)
      return f1 < f2;

    return probs[i - 1] < probs[j - 1];
  });

  return res;
}

// [[Rcpp::export]]
std::vector<bool>
FilterPredictions(const s_vec_t &not_filtered_umis, const s_vec_t &base_umis, const s_vec_t &target_umis) {
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

// TODO: remove export
//' @export
// [[Rcpp::export]]
DataFrame
PrepareClassifierData(const List &reads_per_umi, const List &neighborhood, const NumericVector &umi_probabilities) {
  auto rpu_data_map = parseList(reads_per_umi);
  ClassifierData res_data(true, parseVector(umi_probabilities));
  auto const &umis = as<s_vec_t>(neighborhood.names());

  for (int i = 0; i < neighborhood.size(); ++i) {
    const auto &umi1_info = UmiInfo(umis[i], rpu_data_map.at(umis[i]));
    for (auto const &umi2 : as<s_vec_t>(as<StringVector>(neighborhood[i]))) {
      res_data.add_umis(umi1_info, UmiInfo(umi2, rpu_data_map.at(umi2)));
    }
  }

  return res_data.to_data_frame();
}

// TODO: remove export
//' @export
// [[Rcpp::export]]
DataFrame PrepareClassifierTrainingData(const List &reads_per_umi_pairs) {
  ClassifierData res_data(false);
  for (const List &rpu_pair : reads_per_umi_pairs) {
    auto umis = as<s_vec_t>(as<StringVector>(rpu_pair.names()));
    res_data.add_umis(UmiInfo(umis[0], as<List>(rpu_pair[0])), UmiInfo(umis[1], as<List>(rpu_pair[1])));
  }

  return res_data.to_data_frame();
}
