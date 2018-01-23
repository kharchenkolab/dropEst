#include <dropestr.h>
#include <umi_processing.h>
#include <ClassifierData.h>

using namespace Rcpp;

ClassifierData::ClassifierRow ClassifierData::fill_row(const UmiInfo &umi_small, const UmiInfo &umi_large) {
  if (umi_small.umi().length() != umi_large.umi().length())
    stop("UMIs must have the same length");

  int diff_pos = -1;
  for (int i = 0; i < umi_small.umi().length(); ++i) {
    char n1 = umi_small.umi()[i], n2 = umi_large.umi()[i];
    if (n1 == n2)
      continue;

    diff_pos = i;
    break;
  }

  if (diff_pos == -1) {
    stop("Equal UMIs are compared");
  }

  return ClassifierRow(umi_small.umi(), umi_large.umi(), umi_small.reads_per_umi(), umi_large.reads_per_umi(),
                       umi_small.quality()[diff_pos]);
}

bool ClassifierData::add_umis(const UmiInfo &umi1, const UmiInfo &umi2) {
  const UmiInfo &umi_small = umi1.reads_per_umi() <= umi2.reads_per_umi() ? umi1 : umi2;
  const UmiInfo &umi_large = umi1.reads_per_umi() <= umi2.reads_per_umi() ? umi2 : umi1;

  if (!this->_presented_umi_pairs.insert(std::make_pair(umi_small.umi(), umi_large.umi())).second)
    return false;

  this->_data.push_back(this->fill_row(umi_small, umi_large));

  return true;
}

DataFrame ClassifierData::to_data_frame() const {
  if (this->_data.empty())
    return DataFrame();

  StringVector umi_small(this->_data.size());
  StringVector umi_large(this->_data.size());
  IntegerVector reads_per_umi_small(this->_data.size());
  IntegerVector reads_per_umi_large(this->_data.size());
  NumericVector mean_log_quality_small(this->_data.size());

  for (size_t i = 0; i < this->_data.size(); ++i) {
    auto const &umi_info = this->_data[i];
    umi_small[i] = umi_info.umi_small;
    umi_large[i] = umi_info.umi_large;
    reads_per_umi_small[i] = umi_info.reads_per_umi_small;
    reads_per_umi_large[i] = umi_info.reads_per_umi_large;
    mean_log_quality_small[i] = umi_info.log_quality_small;
  }

  return DataFrame::create(_["Base"] = umi_small, _["Target"] = umi_large, _["MinRpU"] = reads_per_umi_small,
                           _["MaxRpU"] = reads_per_umi_large, _["Quality"] = mean_log_quality_small);
}

ClassifierData::ClassifierRow::ClassifierRow(const std::string &umi_small, const std::string umi_large,
                                             int reads_per_umi_small, int reads_per_umi_large, double log_quality_small)
  : umi_small(umi_small)
  , umi_large(umi_large)
  , reads_per_umi_small(reads_per_umi_small)
  , reads_per_umi_large(reads_per_umi_large)
  , log_quality_small(log_quality_small)
{}
