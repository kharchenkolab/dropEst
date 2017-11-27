#include <dropestr.h>
#include <umi_processing.h>
#include <ClassifierData.h>
#include <CppMap.h>

#include <cmath>

using namespace Rcpp;

ClassifierData::ClassifierData(const sd_map_t &umi_probabilities, const sd_map_t &probability_normalizers)
  : _umi_probabilities_map(umi_probabilities)
  , _probability_normalizers(probability_normalizers)
  , _save_probs(!_umi_probabilities_map.empty())
{}

ClassifierData::ClassifierRow ClassifierData::fill_row(const UmiInfo &umi_small, const UmiInfo &umi_large, double umi_prob) {
  if (umi_small.umi().length() != umi_large.umi().length())
    stop("UMIs must have the same length");

  int diff_pos = -1;
  std::string nuc_diff = "NN";
  for (int i = 0; i < umi_small.umi().length(); ++i) {
    char n1 = umi_small.umi()[i], n2 = umi_large.umi()[i];
    if (n1 == n2)
      continue;

    diff_pos = i;
    nuc_diff = std::string({n1, n2});
    break;
  }

  if (diff_pos == -1) {
    stop("Equal UMIs are compared");
  }

  return ClassifierRow(umi_small.umi(), umi_large.umi(), umi_small.reads_per_umi(), umi_large.reads_per_umi(),
                       umi_large.umi()[diff_pos], NUCL_PAIR_INDS.at(nuc_diff), diff_pos, umi_prob,
                       umi_small.quality()[diff_pos]);
}

bool ClassifierData::add_umis(const UmiInfo &umi1, const UmiInfo &umi2) {
  const UmiInfo &umi_small = umi1.reads_per_umi() <= umi2.reads_per_umi() ? umi1 : umi2;
  const UmiInfo &umi_large = umi1.reads_per_umi() <= umi2.reads_per_umi() ? umi2 : umi1;

  if (!this->_presented_umi_pairs.insert(std::make_pair(umi_small.umi(), umi_large.umi())).second)
    return false;

  double umi_prob = -1;
  if (this->_save_probs) {
    umi_prob = this->_umi_probabilities_map.at(umi_small.umi());
    if (!this->_probability_normalizers.empty()) {
      umi_prob /= this->_probability_normalizers.at(umi_large.umi());
    }
  }
  this->_data.push_back(this->fill_row(umi_small, umi_large, umi_prob));
  return true;
}

DataFrame ClassifierData::to_data_frame() const {
  if (this->_data.empty())
    return DataFrame();

  StringVector umi_small(this->_data.size());
  StringVector umi_large(this->_data.size());
  StringVector nucleotide_large(this->_data.size());
  IntegerVector reads_per_umi_small(this->_data.size());
  IntegerVector reads_per_umi_large(this->_data.size());
  IntegerVector nucleotide_substitution(this->_data.size());
  IntegerVector difference_position(this->_data.size());
  NumericVector umi_prob(this->_data.size());
  NumericVector mean_log_quality_small(this->_data.size());

  for (size_t i = 0; i < this->_data.size(); ++i) {
    auto const &umi_info = this->_data[i];
    umi_small[i] = umi_info.umi_small;
    umi_large[i] = umi_info.umi_large;
    nucleotide_large[i] = std::string(1, umi_info.nucleotide_large);
    reads_per_umi_small[i] = umi_info.reads_per_umi_small;
    reads_per_umi_large[i] = umi_info.reads_per_umi_large;
    nucleotide_substitution[i] = umi_info.nucleotide_substitution;
    difference_position[i] = umi_info.difference_position;
    umi_prob[i] = std::log(umi_info.umi_prob);
    mean_log_quality_small[i] = umi_info.log_quality_small;
  }

  auto res = DataFrame::create(_["Base"] = umi_small, _["Target"] = umi_large, _["NucleotideLarge"] = nucleotide_large,
                               _["Position"] = difference_position, _["Nucleotides"] = nucleotide_substitution,
                               _["MinRpU"] = reads_per_umi_small, _["MaxRpU"] = reads_per_umi_large,
                               _["Quality"] = mean_log_quality_small);

  if (this->_save_probs) {
    res["UmiProb"] = umi_prob;
  }

  return res;
}

ClassifierData::ClassifierRow::ClassifierRow(const std::string &umi_small, const std::string umi_large,
                                             int reads_per_umi_small, int reads_per_umi_large, char nucleotide_large,
                                             int nucleotide_substitution, int difference_position, double umi_prob,
                                             double log_quality_small)
  : umi_small(umi_small)
  , umi_large(umi_large)
  , reads_per_umi_small(reads_per_umi_small)
  , reads_per_umi_large(reads_per_umi_large)
  , nucleotide_large(nucleotide_large)
  , nucleotide_substitution(nucleotide_substitution)
  , difference_position(difference_position)
  , umi_prob(umi_prob)
  , log_quality_small(log_quality_small)
{}
