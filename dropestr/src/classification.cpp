#include <dropestr.h>
#include <umi_processing.h>
#include <ClassifierData.h>

using namespace Rcpp;

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

//// [[Rcpp::export]]
//NumericVector PredictLeftPart(const List &classifier, const DataFrame &predict_data, int gene_size) {
//  using NV = NumericVector;
//  using IV = IntegerVector;
//
//  List error_distributions = as<List>(classifier["Negative"]);
//  List common_distributions = as<List>(classifier["Common"]);
//
//  NV nucl_prob_err = log(as<NV>(as<NV>(error_distributions["Nucleotides"])[as<IV>(predict_data["Nucleotides"])]));
//  NV position_prob_err = log(as<NV>(as<NV>(error_distributions["Position"])[as<IV>(predict_data["Position"])]));
//
//  NV min_rpu_prob_err = std::log(as<double>(error_distributions["MinRpU"])) * as<NV>(predict_data["MinRpU"]);
//  NV max_rpu_prob_err = std::log(as<double>(error_distributions["MaxRpU"])) * as<NV>(predict_data["MaxRpU"]);
//
//  IV quantized_quality = Quantize(as<NV>(predict_data["Quality"]), as<NV>(classifier["QualityQuantBorders"]));
//
//  NV rpu_prob = as<NV>(common_distributions["RpuProbs"]);
//
//  NV min_rpu_prob(max_rpu_prob_err.size());
//  NV max_rpu_prob(max_rpu_prob_err.size());
//  NV quality_prob(max_rpu_prob_err.size());
//  NV quality_prob_err(max_rpu_prob_err.size());
//
//  for (int i = 0; i < max_rpu_prob_err.size(); ++i) {
//    min_rpu_prob[i] = std::log(rpu_prob.at(as<IV>(predict_data["MinRpU"])[i] - 1));
//    max_rpu_prob[i] = std::log(rpu_prob.at(as<IV>(predict_data["MaxRpU"])[i] - 1));
//
//    int current_quality_quant = quantized_quality[i];
//    quality_prob[i] = std::log(as<NV>(common_distributions["Quality"]).at(current_quality_quant));
//    quality_prob_err[i] = std::log(as<NV>(error_distributions["Quality"]).at(current_quality_quant));
//  }
//
//  NV umi_prob_pos = log(NV(1 - vpow(1 - as<NV>(predict_data["UmiProb"]), gene_size)));
//
//  return exp((nucl_prob_err + position_prob_err + min_rpu_prob_err + max_rpu_prob_err) - //  + quality_prob_err
//             (umi_prob_pos + min_rpu_prob + max_rpu_prob)); //  + quality_prob
//}

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

// [[Rcpp::export]]
DataFrame PrepareClassifierData(const List &reads_per_umi) {
    auto rpu_data_map = parseList(reads_per_umi);
    ClassifierData res_data;

    try {
      auto neighbourhood = SubsetAdjacentUmis(as<s_vec_t>(reads_per_umi.names())); // TODO: now it accepts only paired UMIs
      for (auto const &it : neighbourhood) {
        const auto &umi1_info = UmiInfo(it.first, rpu_data_map.at(it.first));
        for (auto const &umi2 : it.second) {
          res_data.add_umis(umi1_info, UmiInfo(umi2, rpu_data_map.at(umi2)));
        }
      }
    }
    catch (std::exception &ex) {
      stop("Can't prepare classifier data: " + std::string(ex.what()));
    }

    return res_data.to_data_frame();
  }

// [[Rcpp::export]]
DataFrame PrepareClassifierTrainingData(const List &reads_per_umi_pairs) {
  ClassifierData res_data;
  for (const List &rpu_pair : reads_per_umi_pairs) {
    auto umis = as<s_vec_t>(as<StringVector>(rpu_pair.names()));
    res_data.add_umis(UmiInfo(umis[0], as<List>(rpu_pair[0])), UmiInfo(umis[1], as<List>(rpu_pair[1])));
  }

  return res_data.to_data_frame();
}
