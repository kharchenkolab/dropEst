#include "dropestr.h"
#include "umi_processing.h"

using namespace Rcpp;

s_vec_t GetUmisList(unsigned umi_len);

//' @export
// [[Rcpp::export]]
SEXP BuildCountMatrix(const List &reads_per_umi_per_cell) {
  Progress p(0, false);
  using Triplet=Eigen::Triplet<unsigned>;

  auto umis_per_gene = as<IntegerVector>(reads_per_umi_per_cell["umis_per_gene"]);
  auto cells = as<CharacterVector>(reads_per_umi_per_cell["cells"]);
  auto genes = as<CharacterVector>(reads_per_umi_per_cell["genes"]);
  auto cell_indexes = as<IntegerVector>(reads_per_umi_per_cell["cell_indexes"]);
  auto gene_indexes = as<IntegerVector>(reads_per_umi_per_cell["gene_indexes"]);

  std::vector<Triplet> triplet_list(umis_per_gene.size());
  for (int i = 0; i < umis_per_gene.size(); ++i)
  {
    if (p.check_abort())
      return IntegerMatrix();

    triplet_list[i] = Triplet(gene_indexes[i], cell_indexes[i], umis_per_gene[i]);
  }

  Eigen::SparseMatrix<unsigned> mat(genes.size(), cells.size());
  mat.setFromTriplets(triplet_list.begin(), triplet_list.end());

  S4 res(wrap(mat));
  res.slot("Dimnames") = List::create(genes, cells);
  return res;
}

//' @export
// [[Rcpp::export]]
List TrimUmis(const List &rpu_per_cell, int trim_length, bool reverse=false) {
  UmisInfo umis_info(rpu_per_cell);
  UmisInfo umis_info_trimmed;

  for (auto const &umi_iter : umis_info.info()) {
    auto const &umi_info = umi_iter.second;
    size_t trim_start = reverse ? (umi_info.umi().length() - trim_length) : 0;
    UmiInfo::quality_t trimmed_quality(umi_info.quality().begin() + trim_start,
                                       umi_info.quality().begin() + trim_start + trim_length);
    UmiInfo trimmed_umi(umi_info.umi().substr(trim_start, trim_length), umi_info.reads_per_umi(), trimmed_quality);
    umis_info_trimmed.add_umi(trimmed_umi);
  }

  return umis_info_trimmed.to_list();
}

//' Estimate a distribution of observed UMI probabilities.
//'
//' @param umis_per_gene_per_cell list of vectors: number of UMIs per gene per cell (zeros can be omitted).
//' @param umi_length length of UMI.
//' @param smooth smooth term, which is added to each UMI probability in case if some UMIs have only few observations.
//' @return Vector of UMI probabilities.
//' @export
// [[Rcpp::export]]
si_map_t GetUmisDistribution(List umis_per_gene_per_cell, int smooth = 1) {
  si_map_t res;

  for (const List &rpus: umis_per_gene_per_cell) {
    for (auto const &umi : GeneInfo(rpus).umis) {
      res[umi]++;
    }
  }

  if (res.empty()) {
    warning("UMI distribution is empty");
    return res;
  }

  unsigned umi_length = res.begin()->first.length();
  for (auto const &umi : GetUmisList(umi_length)) {
    res[umi] += smooth;
  }

  return res;
}

void addNucleotideToUmi(char *current_umi, unsigned current_pos, unsigned umi_len, s_vec_t &umis) {
  if (current_pos > umi_len) {
    umis.push_back(std::string(current_umi, umi_len));
    return;
  }

  for (char n : NUCLEOTIDES) {
    current_umi[current_pos] = n;
    addNucleotideToUmi(current_umi, current_pos + 1, umi_len, umis);
  }
}

// [[Rcpp::export]]
s_vec_t GetUmisList(unsigned umi_len) {
  s_vec_t res;
  char *umi = new char[umi_len];
  addNucleotideToUmi(umi, 0, umi_len, res);
  delete[] umi;

  return res;
}

// [[Rcpp::export]]
List ConcatLists(const List &lists) {
  int total_len = 0;
  for (auto const &lst : lists) {
    total_len += as<List>(lst).size();
  }

  List res(total_len);
  int res_ind = 0;
  for (auto const &lst : lists) {
    for (auto const &val : as<List>(lst)) {
      res[res_ind++] = val;
    }
  }
  return res;
}

// [[Rcpp::export]]
std::unordered_set<int> GetMirrorPairs(const StringVector &pairs, const NumericVector &probs, double tol=1e-5) {
  std::unordered_set<int> answers;

  std::map<ss_pair, int> pair_inds;

  for (int i = 0; i < as<IntegerVector>(pairs.attr("dim"))[0]; ++i) {
    auto pair = std::make_pair(as<std::string>(pairs(i, 0)), as<std::string>(pairs(i, 1)));
    pair_inds.emplace(pair, i);
  }

  for (int i = 0; i < as<IntegerVector>(pairs.attr("dim"))[0]; ++i) {
    auto reverse_pair = std::make_pair(as<std::string>(pairs(i, 1)), as<std::string>(pairs(i, 0)));
    auto reverse_ind_it = pair_inds.find(reverse_pair);
    if (reverse_ind_it != pair_inds.end()) {
      double reverse_prob = probs.at(reverse_ind_it->second), direct_prob = probs.at(i);
      if (reverse_prob > direct_prob + tol) {
        answers.insert(i);
      }
      else if (reverse_prob < direct_prob - tol) {
        answers.insert(reverse_ind_it->second);
      }
      else {
        answers.insert(std::min(reverse_ind_it->second, i));
      }
    }
  }

  return answers;
}

// [[Rcpp::export]]
unsigned NumberOfNucleotidePairs() {
  return unsigned(NUCL_PAIR_INDS.size());
}
