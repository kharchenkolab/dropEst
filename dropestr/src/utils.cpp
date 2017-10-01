#include "dropestr.h"
#include "umi_processing.h"

using namespace Rcpp;

s_vec_t GetUmisList(unsigned umi_len);

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

//' Parse UMIgs.
//'
//' @param reads_per_umigs data from the estimation step.
//' @param umi_length length of UMI.
//' @return List of lists of vectors with number of reads per UMI per gene per cell.
//'
//' @export
// [[Rcpp::export]]
List ParseUmisPerGene(const List &reads_per_umigs, int umi_length) {
  List res(reads_per_umigs.size());

  auto const &cell_names = as<s_vec_t>(as<StringVector>(reads_per_umigs.names()));
  for (int cell_ind = 0; cell_ind < reads_per_umigs.size(); ++cell_ind) {
    ssi_map_t cur_genes;
    const IntegerVector &rp_umigs = as<IntegerVector>(reads_per_umigs[cell_ind]);
    const auto &umigs = as<s_vec_t>(as<StringVector>(rp_umigs.names()));

    for (int i = 0; i < rp_umigs.size(); ++i) {
      const auto &umig = umigs[i];
      std::string umi(umig.substr(0, umi_length));
      if (umi.find('N') != std::string::npos)
        continue;

      std::string gene(umig.substr(umi_length));
      cur_genes[gene][umi] = rp_umigs[i];
    }

    res[cell_ind] = cur_genes;
  }

  res.attr("names") = as<StringVector>(reads_per_umigs.names());
  return res;
}

// [[Rcpp::export]]
List TrimUmis(const List &rpu_per_cell, int trim_length) {
  List res(rpu_per_cell.size());

  auto const &cell_names = as<s_vec_t>(as<StringVector>(rpu_per_cell.names()));
  for (int cell_ind = 0; cell_ind < rpu_per_cell.size(); ++cell_ind) {
    const auto &rpu_per_gene = as<List>(rpu_per_cell[cell_ind]);
    List cell_rpus(rpu_per_gene.size());

    for (int gene_ind = 0; gene_ind < rpu_per_gene.size(); ++gene_ind) {
      auto const &reads_per_umi = as<IntegerVector>(rpu_per_gene[gene_ind]);
      auto const &umis = as<s_vec_t>(as<StringVector>(reads_per_umi.names()));
      si_map_t trimmed_rpus;
      for (int umi_ind = 0; umi_ind < reads_per_umi.size(); ++umi_ind) {
        trimmed_rpus[umis[umi_ind].substr(0, trim_length)] += reads_per_umi[umi_ind];
      }

      cell_rpus[gene_ind] = wrap(trimmed_rpus);
    }

    cell_rpus.attr("names") = as<StringVector>(rpu_per_gene.names());
    res[cell_ind] = cell_rpus;
  }

  res.attr("names") = cell_names;
  return res;
}

//' @export
// [[Rcpp::export]]
List AddIndexesToRpU(const List &reads_per_umi_per_cb, const std::vector<std::string> &umis) { // TODO: not export
  si_map_t umi_inds;
  for (int i = 0; i < umis.size(); ++i) {
    umi_inds.emplace(umis[i], i);
  }

  List res(reads_per_umi_per_cb.size());
  for (int rpu_ind = 0; rpu_ind < reads_per_umi_per_cb.size(); ++rpu_ind) {
    auto const &rpus = as<List>(reads_per_umi_per_cb[rpu_ind]);
    auto const cur_umis = as<s_vec_t>(GeneInfo(rpus).umis);
    IntegerVector indexes(cur_umis.size());
    for (int umi_ind = 0; umi_ind < indexes.size(); ++umi_ind) {
      indexes[umi_ind] = umi_inds.at(cur_umis[umi_ind]);
    }

    res[rpu_ind] = List::create(_["rpus"] = rpus, _["indexes"] = indexes);
  }

  return res;
}

//' Estimate a distribution of observed UMI probabilities.
//'
//' @param umis_per_gene_per_cell list of vectors: number of UMIs per gene per cell (zeros can be omitted).
//' @param umi_length length of UMI.
//' @param smooth smooth term, which is added to each UMI probability in case if some UMIs have only few observations.
//' @return Vector of UMI probabilities.
//'
//' @export
// [[Rcpp::export]]
si_map_t GetUmisDistribution(List umis_per_gene_per_cell, int smooth = 1) {
  si_map_t res;

  for (const List &rpus: umis_per_gene_per_cell) {
    for (auto const &umi : as<s_vec_t>(GeneInfo(rpus).umis)) {
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
