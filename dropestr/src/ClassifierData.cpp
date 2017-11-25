#include <dropestr.h>
#include <umi_processing.h>
#include <ClassifierData.h>
#include <CppMap.h>

using namespace Rcpp;

ClassifierData::ClassifierData(bool force_neighbours, const sd_map_t &umi_probabilities_map)
  : _umi_probabilities_map(umi_probabilities_map)
  , _save_probs(!_umi_probabilities_map.empty())
  , _force_neighbours(force_neighbours)
{}

ClassifierData::ClassifierRow ClassifierData::fill_row(const UmiInfo &umi_small, const UmiInfo &umi_large, double umi_prob) {
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

bool ClassifierData::add_umis(const UmiInfo &umi1, const UmiInfo &umi2) {
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

DataFrame ClassifierData::to_data_frame() const {
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

ClassifierData::ClassifierRow::ClassifierRow(const std::string &umi_small, const std::string umi_large,
                                             int reads_per_umi_small, int reads_per_umi_large, int edit_distance,
                                             int nucleotide_substitution, int difference_position, double umi_prob,
                                             double log_quality_small)
  : umi_small(umi_small)
  , umi_large(umi_large)
  , reads_per_umi_small(reads_per_umi_small)
  , reads_per_umi_large(reads_per_umi_large)
  , edit_distance(edit_distance)
  , nucleotide_substitution(nucleotide_substitution)
  , difference_position(difference_position)
  , umi_prob(umi_prob)
  , log_quality_small(log_quality_small)
{}
