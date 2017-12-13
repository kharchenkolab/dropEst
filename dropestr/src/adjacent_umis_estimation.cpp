#include "dropestr.h"
#include "umi_processing.h"

using namespace Rcpp;

template<typename T>
s_vec_t getAdjacentUmis(const std::string &umi, const T &filter, bool need_filter) {
  static_assert((std::is_same<s_set_t, T>::value || std::is_same<si_map_t, T>::value), "Only map or set");
  s_vec_t res;
  for (int i = 0; i < umi.length(); ++i) {
    std::string cur_umi(umi);
    for (char n : NUCLEOTIDES) {
      if (umi[i] == n)
        continue;

      cur_umi[i] = n;

      if (!need_filter || filter.find(cur_umi) != filter.end()) {
        res.push_back(cur_umi);
      }
    }
  }

  return res;
}

//' @export
// [[Rcpp::export]]
s_vec_t GetAdjacentUmis(const std::string &umi) { // TODO: not @export
  return getAdjacentUmis(umi, s_set_t(), false);
}

// [[Rcpp::export]]
std::vector<bool> GetCrossmergedMask(const s_vec_t &base_umis, const s_vec_t &target_umis) {
  si_map_t targets_num = ValueCountsC(target_umis);
  std::vector<bool> res_mask(base_umis.size(), false);
  for (int i = 0; i < base_umis.size(); ++i) {
    auto iter = targets_num.find(base_umis[i]);
    if (iter != targets_num.end() && iter->second == 1) {
      res_mask[i] = true;
    }
  }

  return res_mask;
}

// [[Rcpp::export]]
std::vector<bool> ResolveUmiDependencies(const s_vec_t &base_umis, const s_vec_t &target_umis, bool verbose=false) { // TODO: remove verbose
  if (base_umis.size() != target_umis.size())
    stop("All vectors must have the same size");

  si_map_t inds_by_base;
  s_vec_t base_uniq; // TODO: remove
  for (int i = 0; i < base_umis.size(); ++i) {
    auto iter = inds_by_base.emplace(base_umis[i], inds_by_base.size());
    if (iter.second) {
      base_uniq.push_back(base_umis[i]);
    }
  }

  std::vector<int> merge_targets(inds_by_base.size());
  std::iota(merge_targets.begin(), merge_targets.end(), 0);

  for (size_t id = 0; id < base_umis.size(); ++id) {
    int base_umi_id = inds_by_base.at(base_umis[id]); // if several umis with the same base_umi are presented
    if (merge_targets[base_umi_id] != base_umi_id)
      continue;

    if (verbose) {
      std::cout << id << ", " << base_umi_id << ": " << base_uniq[base_umi_id] << "\n";
    }
    auto const &target_umi = target_umis[id];
    auto iter = inds_by_base.find(target_umi);
    int target_id = (iter == inds_by_base.end()) ? -1 : iter->second;
    if (verbose) {
      std::cout << "\t" << target_umis[id] << " (" << target_id << ")\n";
    }
    while (target_id != -1 && target_id != base_umi_id && target_id != merge_targets[target_id]) {
      target_id = merge_targets[target_id];
      if (verbose) {
        std::cout << "\t" << ((target_id != -1) ? base_uniq[target_id] : std::to_string(-1)) << " (" << target_id << ")\n";
      }
    }

    merge_targets[base_umi_id] = target_id;
  }

  std::vector<bool> is_filtered;
  for (int i = 0; i < base_umis.size(); ++i) {
    int base_umi_id = inds_by_base.at(base_umis[i]);
    is_filtered.push_back(merge_targets[base_umi_id] != base_umi_id);
    if (verbose) {
      std::cout << "(" << base_umi_id << ", " << merge_targets[base_umi_id] << ") ";
    }
  }
  if (verbose) {
    std::cout << std::endl;
  }

  return is_filtered;
}

// [[Rcpp::export]]
std::unordered_map<std::string, s_vec_t> SubsetAdjacentUmis(const s_vec_t &umis) { // TODO: not export
  s_set_t umis_set(umis.begin(), umis.end());
  std::unordered_map<std::string, s_vec_t> res;

  for (auto const &umi : umis) {
    res[umi] = getAdjacentUmis(umi, umis_set, true);
  }

  return res;
}

//' Fill information about adjacent UMIs, their probabilities and differences for each UMI.
//'
//' @param umi_probabilites vector of UMI probabilities.
//' @param adjacent_only logical, return only the list of adjacent UMIs.
//' @param show_progress show progress bar.
//' @return List with the information about adjacent UMIs.
//'
//' @export
// [[Rcpp::export]]
List FillAdjacentUmisData(const NumericVector &umi_probabilites, bool adjacent_only=false, bool show_progress=false) {
  auto const &umis = as<s_vec_t>(umi_probabilites.names());
  const sd_map_t umi_probs_map = parseVector(umi_probabilites);

  sd_map_t adjacent_probs;
  Progress progress(umis.size(), show_progress);
  std::unordered_map<std::string, s_vec_t> neighbours;
  for (auto const &cur_umi : umis) {
    auto cur_neighbours = neighbours.insert(std::make_pair(cur_umi, s_vec_t())).first;
    double sum_prob = 0;
    for (int i = 0; i < cur_umi.length(); ++i) {
      std::string neighb_umi(cur_umi);
      double cur_umi_prob = umi_probs_map.at(cur_umi);
      for (char n : NUCLEOTIDES) {
        if (cur_umi[i] == n)
          continue;

        neighb_umi[i] = n;
        cur_neighbours->second.push_back(neighb_umi);

        if (adjacent_only)
          continue;

        double neighb_umi_prob = umi_probs_map.at(neighb_umi);
        sum_prob += neighb_umi_prob;
      }
    }

    if (!adjacent_only) {
      adjacent_probs.emplace(cur_umi, sum_prob);
    }

    if (progress.check_abort())
      return List();

    progress.increment();
  }

  if (adjacent_only)
    return wrap(neighbours);

  return List::create(_["adjacent.umis"] = neighbours, _["probabilities"] = adjacent_probs);
}

//' @export
// [[Rcpp::export]]
List GetAdjacentUmisNum(const IntegerVector &reads_per_umi_from, const IntegerVector &reads_per_umi_to) {
  IntegerVector total_neighbours_num(reads_per_umi_from.size(), 0),
                smaller_neighbours_num(reads_per_umi_from.size(), 0),
                larger_neighbours_num(reads_per_umi_from.size(), 0);

  auto const &umis = as_s_vec(reads_per_umi_from.names());
  s_set_t umis_set(umis.begin(), umis.end());

  si_map_t rpu_map_to = parseVector(reads_per_umi_to);

  for (size_t i = 0; i < umis.size(); ++i) {
    int cur_rpu = reads_per_umi_from[i];

    int total_nn = 0, smaller_nn = 0, larger_nn = 0;
    for (auto const &neighbour : getAdjacentUmis(umis[i], umis_set, true)) {
      auto neighb_rpu_it = rpu_map_to.find(neighbour);
      if (neighb_rpu_it == rpu_map_to.end())
        continue;

      if (neighb_rpu_it->second <= cur_rpu) {
        smaller_nn++;
      }
      else {
        larger_nn++;
      }

      total_nn++;
    }

    smaller_neighbours_num[i] = smaller_nn;
    larger_neighbours_num[i] = larger_nn;
    total_neighbours_num[i] = total_nn;
  }


  smaller_neighbours_num.attr("names") = wrap(umis);
  larger_neighbours_num.attr("names") = wrap(umis);
  total_neighbours_num.attr("names") = wrap(umis);

  return List::create(_["Smaller"]=smaller_neighbours_num,
                      _["Larger"]=larger_neighbours_num,
                      _["Total"]=total_neighbours_num);
}


//' @export
// [[Rcpp::export]]
NumericMatrix FillDpMatrix(double prior_prob, int neighbours_num, int max_umi_per_cell) { //TODO: not export this
  int n_col = max_umi_per_cell, n_row = neighbours_num + 1;
  NumericMatrix dp_matrix(n_row, n_col);

  double prod = 1;
  for (int i = 0; i < n_col; ++i) {
    dp_matrix(0, i) = prod;
    prod *= 1 - prior_prob;
  }

  for (int row = 1; row < n_row; ++row) {
    for (int col = 1; col < n_col; ++col) {
      dp_matrix(row, col) = dp_matrix.at(row - 1, col - 1) * prior_prob * (1.0 - double(row - 1.0) / neighbours_num) +
        dp_matrix.at(row, col - 1) * (1.0 - prior_prob * (1.0 - double(row) / neighbours_num));
    }
  }

  return dp_matrix;
}

void fillReverseCumSum(int max_neighbour_num, int larger_nn, int umi_ind, const NumericVector &distr, NumericMatrix &res) {
  double cum_sum = 0;
  for (int nn = max_neighbour_num; nn >= larger_nn; --nn) {
    cum_sum += distr[nn];
    res(nn - larger_nn, umi_ind) = cum_sum;
  }
}

void fillCumSumRatio(int max_neighbour_num, int smaller_nn, int larger_nn, int umi_ind, const NumericVector &distr,
                     NumericMatrix &res, bool log_probs=false, double tol=1e-20) {
  if (distr.size() <= larger_nn + smaller_nn) {
    stop("Strange numbers: max=" + std::to_string(max_neighbour_num) + ", larget=" + std::to_string(larger_nn) +
      ", smaller=" + std::to_string(smaller_nn));
  }

  double reverse_cum_sum = 0, direct_cum_sum = 0;
  std::vector<double> direct_cum_sums(max_neighbour_num + 1, 0); // To prevent underflow
  for (int nn = larger_nn; nn <= larger_nn + smaller_nn; ++nn) {
    direct_cum_sum += distr[nn];
    direct_cum_sums[nn] = direct_cum_sum;
  }

  if (direct_cum_sum < tol) {
    for (int nn = larger_nn; nn <= larger_nn + smaller_nn; ++nn) {
      direct_cum_sums[nn] = 1e10; // Any huge number
    }
  }

  for (int nn = larger_nn + smaller_nn; nn > larger_nn; --nn) {
    reverse_cum_sum += distr[nn];
    double log_value = std::log(reverse_cum_sum) - std::log(direct_cum_sums[nn - 1]);
    res(nn - larger_nn, umi_ind) = log_probs ? log_value : std::exp(log_value);
  }
}

// [[Rcpp::export]]
NumericMatrix GetSmallerNeighboursDistributionsBySizes(const List &dp_matrices, const IntegerVector &larger_neighbours_num,
                                                       const s_vec_t &neighbour_prob_inds, int size_adj, int max_neighbour_num,
                                                       const IntegerVector &smaller_neighbours_num = IntegerVector(), bool log_probs=false,
                                                       bool return_raw=false) {
  if (size_adj == 0)
    stop("Zero gene size");

  // reads_per_umi, larger_neighbours_num and neighbour_prob_inds must have the same order
  si_map_t dp_matrices_index;
  for (int i = 0; i < dp_matrices.size(); ++i) {
    dp_matrices_index.emplace(as<s_vec_t>(dp_matrices.names())[i], i);
  }

  Progress p(0, false);

  NumericMatrix res(max_neighbour_num + 1, larger_neighbours_num.size());
  for (int umi_ind = 0; umi_ind < larger_neighbours_num.size(); ++umi_ind) {
    if (p.check_abort())
      return res;

    auto const &matrix = as<NumericMatrix>(dp_matrices.at(dp_matrices_index.at(neighbour_prob_inds.at(umi_ind))));
    auto const &distr = matrix.column(size_adj - 1);

    int larger_nn = larger_neighbours_num.at(umi_ind);
    double prob_sum = 0;
    if (distr.size() <= max_neighbour_num)
      stop("Too small distr size");

    for (int nn = larger_nn; nn <= max_neighbour_num; ++nn) {
      prob_sum += distr[nn];
    }

    if (prob_sum < 1e-10) {
      prob_sum = 1;
    }

    if (return_raw) {
      for (int nn = larger_nn; nn <= max_neighbour_num; ++nn) {
        res(nn - larger_nn, umi_ind) = distr[nn] / prob_sum;
      }

      continue;
    }

    if (smaller_neighbours_num.size() != 0) {
      fillCumSumRatio(max_neighbour_num, smaller_neighbours_num.at(umi_ind), larger_nn, umi_ind, distr / prob_sum, res, log_probs);
      continue;
    }

    fillReverseCumSum(max_neighbour_num, larger_nn, umi_ind, distr / prob_sum, res);
  }

  colnames(res) = as<StringVector>(larger_neighbours_num.names());

  return res;
}

// [[Rcpp::export]]
List FilterUmisInGeneClassic(const List &reads_per_umi, double mult=1) {
  if (reads_per_umi.size() == 1)
    return reads_per_umi;

  Progress p(0, false);
  const GeneInfo info(reads_per_umi);
  si_map_t rpu_map;
  for (int i = 0; i < info.umis.size(); ++i) {
    rpu_map.emplace(info.umis[i], info.reads_per_umi[i]);
  }

  auto const neighbourhood = SubsetAdjacentUmis(info.umis);
  si_map_t umi_inds;
  for (int i = 0; i < info.umis.size(); ++i) {
    umi_inds[info.umis[i]] = i;
  }

  std::vector<std::string> base_umis, target_umis;
  for (int i = 0; i < info.umis.size(); ++i) {
    auto const &neighbours = neighbourhood.at(info.umis[i]);
    int cur_rpu = rpu_map.at(info.umis[i]);

    for (auto const &neighbour : neighbours) {
      if (rpu_map.at(neighbour) < cur_rpu * mult - EPS) // rpu_neighb >= cur_rpu
        continue;

      base_umis.push_back(info.umis[i]);
      target_umis.push_back(neighbour);
    }

    if (p.check_abort())
      break;
  }

  std::vector<bool> is_base_filt = ResolveUmiDependencies(base_umis, target_umis);
  LogicalVector filt_mask(reads_per_umi.size(), true);
  for (int i = 0; i < is_base_filt.size(); ++i) {
    if (is_base_filt[i]) {
      filt_mask[umi_inds[base_umis[i]]] = false;
    }
  }

  return reads_per_umi[filt_mask];
}
