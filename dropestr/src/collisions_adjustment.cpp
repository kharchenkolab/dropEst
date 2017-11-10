#include "dropestr.h"
#include "umi_processing.h"

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector GetTrimCollisionsNum(const List &rpu_per_gene, int trim_length) {
  IntegerVector res(rpu_per_gene.size());

  for (int gene_ind = 0; gene_ind < rpu_per_gene.size(); ++gene_ind) {
    const GeneInfo gene_info(as<List>(rpu_per_gene[gene_ind]));
    s_set_t trimmed_rpus;

    int collisions_num = 0;
    for (int umi_ind = 0; umi_ind < gene_info.reads_per_umi.size(); ++umi_ind) {
      auto ins = trimmed_rpus.insert(gene_info.umis[umi_ind].substr(0, trim_length));
      if (!ins.second) {
        collisions_num++;
      }
    }

    res[gene_ind] = collisions_num;
  }

  return res;
}

class CollisionsAdjuster {
public:
  typedef std::vector<unsigned> size_vec_t;
  typedef std::vector<double> probs_vec_t;

private:
  size_vec_t _adjusted_sizes;
  probs_vec_t _umi_probabilities;
  probs_vec_t _umi_probabilities_neg_prod;

  double _sum_collisions;
  size_t _last_total_gene_size;

private:
  double fpow(double base, long exp) const {
    if (exp == 1)
      return base;

    double result = 1;
    while (exp) {
      if (exp & 1) {
        result *= base;
      }
      exp >>= 1;
      base *= base;
    }

    return result;
  }

  void update_adjusted_sizes(size_t max_gene_expression) {
    for (size_t s = this->_adjusted_sizes.size() + 1; s <= max_gene_expression; ++s) {
      const size_t total_size = s + size_t(this->_sum_collisions);
      double new_umi_prob = 0;
      for (int i = 0; i < this->_umi_probabilities.size(); ++i) {
        this->_umi_probabilities_neg_prod[i] *= this->fpow(1 - this->_umi_probabilities[i], total_size - this->_last_total_gene_size);
        new_umi_prob += this->_umi_probabilities[i] * (1 - this->_umi_probabilities_neg_prod[i]);
      }

      this->_last_total_gene_size = total_size;

      const double collision_num = 1.0 / (1.0 - new_umi_prob) - 1.0;
      this->_sum_collisions += collision_num;
      this->_adjusted_sizes.push_back(std::lround(s + this->_sum_collisions));
    }
  }

public:
  CollisionsAdjuster()
    : _sum_collisions(0)
    , _last_total_gene_size(0)
  {}

  void init(const probs_vec_t &umi_probabilities, size_t max_gene_expression = 0) {
    this->_sum_collisions = 0;
    this->_last_total_gene_size = 0;
    this->_umi_probabilities = umi_probabilities;
    this->_umi_probabilities_neg_prod = probs_vec_t(umi_probabilities.size(), 1);
    CollisionsAdjuster::update_adjusted_sizes(max_gene_expression);
  }

  size_vec_t adjusted_sizes() const {
    return this->_adjusted_sizes;
  }
};

// [[Rcpp::export]]
std::vector<unsigned> FillCollisionsAdjustmentInfo(const std::vector<double> &umi_probabilities, unsigned max_umi_per_gene) {
  CollisionsAdjuster adjuster;
  adjuster.init(umi_probabilities, max_umi_per_gene);
  return adjuster.adjusted_sizes();
}

// '@export
// [[Rcpp::export]]
int AdjustGeneExpressionClassic(int value, int umis_number) {
  if (value == umis_number) {
    return 2 * AdjustGeneExpressionClassic(value - 1, umis_number) - AdjustGeneExpressionClassic(value - 2, umis_number);
  }
  return int(std::round(-std::log(1 - value / double(umis_number)) * umis_number));
}

//' Adjust gene expression value for collisions.
//'
//' @param value gene expression value.
//' @param adjusted_sizes vector of adjusted gene sizes for *observed_sizes*.
//' @param observed_sizes vector of quantized ordered obseved gene sizes.
//' @return Adjusted gene expression value.
//'
//' @export
// [[Rcpp::export]]
unsigned AdjustGeneExpression(unsigned value, const std::vector<int> &adjusted_sizes) {
  if (value > adjusted_sizes.size())
    stop("Too large value of gene expression: " + std::to_string(value) + ", max acceptable is " + std::to_string(adjusted_sizes.size()));

  return adjusted_sizes[value - 1];
}
