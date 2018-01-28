#include "dropestr.h"

using namespace Rcpp;

s_vec_t GetUmisList(unsigned umi_len);

si_map_t parseVector(const IntegerVector &vec) {
  si_map_t res;
  StringVector names = as<StringVector>(vec.names());
  for (int i = 0; i < vec.size(); ++i) {
    res.emplace(std::string(names[i]), vec[i]);
  }

  return res;
}

sd_map_t parseVector(const NumericVector &vec) {
  sd_map_t res;
  StringVector names = as<StringVector>(vec.names());
  for (int i = 0; i < vec.size(); ++i) {
    res.emplace(std::string(names[i]), vec[i]);
  }

  return res;
}

NumericVector vpow(double base, const NumericVector& exp) {
  NumericVector res(exp.size());
  for (int i = 0; i < res.size(); ++i) {
    res[i] = std::pow(base, exp[i]);
  }
  return res;
}

NumericVector vpow(const NumericVector& base, double exp) {
  NumericVector res(base.size());
  for (int i = 0; i < res.size(); ++i) {
    res[i] = std::pow(base[i], exp);
  }
  return res;
}

//' @export
// [[Rcpp::export]]
si_map_t ValueCountsC(const s_vec_t &values) {
  si_map_t res;
  for (auto const &value : values) {
    res[value]++;
  }

  return res;
}

//' @export
// [[Rcpp::export]]
std::unordered_map<int, int> ValueCounts(const std::vector<int> &values) {
  std::unordered_map<int, int> res;
  for (int value : values) {
    res[value]++;
  }

  return res;
}

// [[Rcpp::export]]
List GetUmisDifference(const std::string &umi1, const std::string &umi2, int rpu1, int rpu2, bool force_neighbours, double umi_prob) {
  if (umi1.length() != umi2.length())
    stop("UMIs must have the same length");

  int diff_pos = -1, diff_num = 0;
  std::string nuc_diff = "NN";
  for (int i = 0; i < umi1.length(); ++i) {
    char n1 = umi1[i], n2 = umi2[i];
    if (n1 == n2)
      continue;

    diff_num++;
    if (diff_num != 1)
      continue;

    diff_pos = i;
    nuc_diff[0] = std::min(n1, n2);
    nuc_diff[1] = std::max(n1, n2);

    if (force_neighbours)
      break;
  }

  if (diff_num > 1) {
    diff_pos = runif(1, 0, umi1.length())[0];
    char n1, n2;
    do {
      n1 = NUCLEOTIDES[int(runif(1, 0, NUCLEOTIDES_NUM)[0])];
      n2 = NUCLEOTIDES[int(runif(1, 0, NUCLEOTIDES_NUM)[0])];
    }
    while (n1 == n2);

    nuc_diff[0] = std::min(n1, n2);
    nuc_diff[1] = std::max(n1, n2);
  }

  List res = List::create(_["Position"]=diff_pos, _["Nucleotides"]=NUCL_PAIR_INDS.at(nuc_diff), _["ED"]=diff_num,
                       _["MinRpU"]=std::min(rpu1, rpu2), _["MaxRpU"]=std::max(rpu1, rpu2));

  if (umi_prob > 0) {
    res["UmiProb"] = umi_prob;
  }

  return res;
}

// [[Rcpp::export]]
SEXP BuildCountMatrix(const List &umis_per_gene) {
  si_map_t gene_inds;
  StringVector gene_names;

  Progress p(0, false);

  for (auto const &cell : umis_per_gene) {
    const s_vec_t &genes = as<s_vec_t>(as<StringVector>(as<IntegerVector>(cell).names()));
    for (auto const &gene : genes) {
      auto iter = gene_inds.emplace(gene, gene_names.size());
      if (iter.second) {
        gene_names.push_back(gene);
      }

      if (p.check_abort())
        return IntegerMatrix();
    }
  }

  using Triplet=Eigen::Triplet<unsigned>;
  std::vector<Triplet> triplet_list;
  for (int cell_ind = 0; cell_ind < umis_per_gene.size(); ++cell_ind) {
    auto const &upg = as<IntegerVector>(umis_per_gene[cell_ind]);
    const s_vec_t &genes = as<s_vec_t>(as<StringVector>(upg.names()));

    for (int gene_ind = 0; gene_ind < genes.size(); ++gene_ind) {
      triplet_list.push_back(Triplet(gene_inds.at(genes[gene_ind]), cell_ind, upg[gene_ind]));
    }

    if (p.check_abort())
      return IntegerMatrix();
  }

  Eigen::SparseMatrix<unsigned> mat(gene_names.size(), umis_per_gene.size());
  mat.setFromTriplets(triplet_list.begin(), triplet_list.end());

  S4 res(wrap(mat));
  res.slot("Dimnames") = List::create(gene_names, as<StringVector>(umis_per_gene.names()));
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
  for (int cell_ind = 0; cell_ind < reads_per_umi_per_cb.size(); ++cell_ind) {
    auto const &cell_rpu = as<List>(reads_per_umi_per_cb[cell_ind]);
    List cell_res(cell_rpu.size());

    for (int gene_ind = 0; gene_ind < cell_rpu.size(); ++gene_ind) {
      auto const &cur_umis = as<s_vec_t>(as<StringVector>(as<IntegerVector>(cell_rpu[gene_ind]).names()));
      IntegerVector indexes(cur_umis.size());
      for (int umi_ind = 0; umi_ind < indexes.size(); ++umi_ind) {
        indexes[umi_ind] = umi_inds.at(cur_umis[umi_ind]);
      }

      cell_res[gene_ind] = List::create(_["rpus"] = cell_rpu[gene_ind], _["indexes"] = indexes);
    }

    cell_res.attr("names") = as<StringVector>(cell_rpu.names());
    res[cell_ind] = cell_res;
  }

  res.attr("names") = as<StringVector>(reads_per_umi_per_cb.names());
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

  for (const List &umis_per_gene : umis_per_gene_per_cell) {
    for (const IntegerVector &rpus : umis_per_gene) {
      for (const String &umi : as<StringVector>(rpus.names())) {
        res[umi]++;
      }
    }
  }

  if (res.empty()) {
    warning("UMIs distribution is empty");
    return res;
  }

  int umi_length = res.begin()->first.length();
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
