#pragma once

#include "CommonFuncs.hpp"

s_vec_t GetUmisList(unsigned umi_len);

// [[Rcpp::export]]
IntegerMatrix BuildCountMatrix(const List &umis_per_gene) {
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

  IntegerMatrix count_matrix(Dimension(gene_names.size(), umis_per_gene.size()));
  for (int cell_ind = 0; cell_ind < umis_per_gene.size(); ++cell_ind) {
    auto const &upg = as<IntegerVector>(umis_per_gene[cell_ind]);
    const s_vec_t &genes = as<s_vec_t>(as<StringVector>(upg.names()));
    auto cell_column = count_matrix.column(cell_ind);

    for (int gene_ind = 0; gene_ind < genes.size(); ++gene_ind) {
      cell_column[gene_inds.at(genes[gene_ind])] = upg[gene_ind];
    }

    if (p.check_abort())
      return count_matrix;
  }

  colnames(count_matrix) = as<StringVector>(umis_per_gene.names());
  rownames(count_matrix) = gene_names;

  return count_matrix;
}

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

// [[Rcpp::export]]
List AddIndexesToRpU(const List &reads_per_umi_per_cb, const std::vector<std::string> &umis) {
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

// [[Rcpp::export]]
si_map_t GetUmisDistribution(List umis_per_cell, int umi_length = -1, int smooth = 1) {
  si_map_t res;

  for (const List &umis_per_gene : umis_per_cell) {
    for (const IntegerVector &rpus : umis_per_gene) {
      for (const String &umi : as<StringVector>(rpus.names())) {
        res[umi]++;
      }
    }
  }

  if (umi_length > 0) {
    for (auto const &umi : GetUmisList(umi_length)) {
      res[umi] += smooth;
    }
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
//
// // [[Rcpp::export]]
// IntegerMatrix GetNucleotideFrequencies(IntegerVector umi_distribution) {
//   std::unordered_map<char, int> nucl_inds;
//   nucl_inds['A'] = 0;
//   nucl_inds['C'] = 1;
//   nucl_inds['G'] = 2;
//   nucl_inds['T'] = 3;
//   auto umis = as<StringVector>(umi_distribution.names());
//   IntegerMatrix res_mat(NUCLEOTIDES_NUM, std::string(umis[0]).length());
//   for (int i = 0; i < umis.size(); ++i) {
//     int c_ind = 0;
//     int cur_freq = umi_distribution[i];
//     for (char c : std::string(umis[i])) {
//       res_mat(nucl_inds[c], c_ind) += cur_freq;
//       c_ind++;
//     }
//   }
//
//   return res_mat;
// }

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