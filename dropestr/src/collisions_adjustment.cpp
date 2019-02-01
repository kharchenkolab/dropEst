#include "dropestr.h"
#include "umi_processing.h"

using namespace Rcpp;

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

  void update_adjusted_sizes(size_t max_gene_expression, bool verbose) {
    Progress progress(max_gene_expression - this->_adjusted_sizes.size() - 1, verbose);

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

      if (progress.check_abort())
        break;

      progress.increment();
    }
  }

public:
  CollisionsAdjuster()
    : _adjusted_sizes()
  , _umi_probabilities()
  , _umi_probabilities_neg_prod()
  , _sum_collisions(0)
  , _last_total_gene_size(0)
  {}

  void init(const probs_vec_t &umi_probabilities, size_t max_gene_expression = 0, bool verbose=false) {
    this->_sum_collisions = 0;
    this->_last_total_gene_size = 0;
    this->_umi_probabilities = umi_probabilities;
    this->_umi_probabilities_neg_prod = probs_vec_t(umi_probabilities.size(), 1);
    this->update_adjusted_sizes(max_gene_expression, verbose);
  }

  size_vec_t adjusted_sizes() const {
    return this->_adjusted_sizes;
  }
};

//' @export
// [[Rcpp::export]]
std::vector<unsigned> FillCollisionsAdjustmentInfo(const std::vector<double> &umi_probabilities, unsigned max_umi_per_gene, bool verbose=false) {
  CollisionsAdjuster adjuster;
  adjuster.init(umi_probabilities, max_umi_per_gene, verbose);
  return adjuster.adjusted_sizes();
}

//' @export
// [[Rcpp::export]]
int AdjustGeneExpressionUniform(int value, int umis_number) {
  if (value == umis_number) {
    return 2 * AdjustGeneExpressionUniform(value - 1, umis_number) - AdjustGeneExpressionUniform(value - 2, umis_number);
  }
  return int(std::round(-std::log(1 - value / double(umis_number)) * umis_number));
}

// [[Rcpp::export]]
unsigned DeadjustGeneExpression(double gene_expression, const std::vector<unsigned> &adjusted_expressions) {
  if (adjusted_expressions.empty())
    stop("Empty collisions info");

  if (gene_expression < 0)
    stop("Negative gene expression value: " + std::to_string(gene_expression));

  const double EPS = 1e-3;
  if (gene_expression < 1 + EPS)
    return 1;

  auto upper_iter = std::lower_bound(adjusted_expressions.begin(), adjusted_expressions.end(), gene_expression - EPS);
  if (upper_iter == adjusted_expressions.end())
    stop("Too large value of gene expression: " + std::to_string(gene_expression));

  size_t upper_ind = upper_iter - adjusted_expressions.begin();
  if (upper_ind == 0)
    return 1;

  double upper_expr = adjusted_expressions[upper_ind];
  double lower_expr = adjusted_expressions[upper_ind - 1];

  return std::lround(upper_ind + 1 - (upper_expr - gene_expression) / (upper_expr - lower_expr));
}
