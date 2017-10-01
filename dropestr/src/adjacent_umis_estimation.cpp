#include "dropestr.h"

using namespace Rcpp;

template<typename T>
s_vec_t getAdjacentUmis(const std::string &umi, const T &filter) {
  static_assert((std::is_same<s_set_t, T>::value || std::is_same<si_map_t, T>::value), "Only map or set");
  s_vec_t res;
  for (int i = 0; i < umi.length(); ++i) {
    std::string cur_umi(umi);
    for (char n : NUCLEOTIDES) {
      if (umi[i] == n)
        continue;

      cur_umi[i] = n;

      if (filter.find(cur_umi) != filter.end()) {
        res.push_back(cur_umi);
      }
    }
  }

  return res;
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
List SubsetAdjacentUmis(const s_vec_t &umis) {
  s_set_t umis_set(umis.begin(), umis.end());
  List res(umis.size());

  for (int i = 0; i < umis.size(); ++i) {
    res[i] = getAdjacentUmis(umis[i], umis_set);
  }

  res.attr("names") = wrap(umis);
  return res;
}

std::pair<int, int> GetUmisDifference(const std::string &umi1, const std::string &umi2) {
  if (umi1.length() != umi2.length())
    stop("UMIs must have the same length");

  int diff_pos = -1;
  std::string nuc_diff = "NN";
  for (int i = 0; i < umi1.length(); ++i) {
    char n1 = umi1[i], n2 = umi2[i];
    if (n1 == n2)
      continue;

    diff_pos = i;
    nuc_diff[0] = std::min(n1, n2);
    nuc_diff[1] = std::max(n1, n2);

    break;
  }

  return std::make_pair(diff_pos, NUCL_PAIR_INDS.at(nuc_diff));
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
  NumericVector position_probs(umis.at(0).length(), 0), nucl_probs(NUCL_PAIR_INDS.size(), 0);

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

        auto umi_difference = GetUmisDifference(cur_umi, neighb_umi);

        double pair_prob = neighb_umi_prob * cur_umi_prob;
        position_probs[umi_difference.first] += pair_prob;
        nucl_probs[umi_difference.second] += pair_prob;
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

  double sum_pair_probs = sum(position_probs);
  position_probs = position_probs / sum_pair_probs;
  nucl_probs = nucl_probs / sum_pair_probs;

  return List::create(_["adjacent.umis"] = neighbours, _["probabilities"] = adjacent_probs,
                      _["nucl.probabilities"] = nucl_probs, _["position.probabilities"] = position_probs);
}

// [[Rcpp::export]]
List GetAdjacentUmisNum(const IntegerVector &reads_per_umi_from, const IntegerVector &reads_per_umi_to,
                          const List &neighbourhood, bool total = true, bool larger = false, bool smaller = false) {
  // if (smaller && larger)
  //   stop("You can't use both smaller and larger");

  si_map_t total_neighbours_num, smaller_neighbours_num, larger_neighbours_num;

  auto const &umis = as<StringVector>(reads_per_umi_from.names());
  si_map_t rpu_map = parseVector(reads_per_umi_to);

  for (auto const &neighbours : as<List>(neighbourhood[umis])) {
    int cur_rpu = 0;
    const std::string &cur_umi = std::string(umis[neighbours.index]);

    if (smaller || larger) {
      cur_rpu = rpu_map.at(cur_umi);
    }

    int total_nn = 0, smaller_nn = 0, larger_nn = 0;
    for (auto const &neighbour : as<StringVector>(neighbours)) {
      auto neighb_rpu_it = rpu_map.find(std::string(neighbour));
      if (neighb_rpu_it == rpu_map.end())
        continue;

      if (smaller && neighb_rpu_it->second <= cur_rpu) {
        smaller_nn++;
      }

      if (larger && neighb_rpu_it->second > cur_rpu) {
        larger_nn++;
      }

      if (total) {
        total_nn++;
      }
    }

    if (smaller) {
      smaller_neighbours_num.emplace(cur_umi, smaller_nn);
    }

    if (larger) {
      larger_neighbours_num.emplace(cur_umi, larger_nn);
    }

    if (total) {
      total_neighbours_num.emplace(cur_umi, total_nn);
    }
  }

  return List::create(_["Smaller"]=smaller_neighbours_num,
                      _["Larger"]=larger_neighbours_num,
                      _["Total"]=total_neighbours_num);
}


//' @export
// [[Rcpp::export]]
NumericMatrix FillDpMatrix(double prior_prob, int neighbours_num, int max_umi_per_cell) { //TODO: not export this
  int n_col = max_umi_per_cell + 1, n_row = neighbours_num + 1;
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

void fillCumSumRatio(int max_neighbour_num, int smaller_nn, int larger_nn, int umi_ind, const NumericVector &distr, NumericMatrix &res, double tol=1e-10) {
  double reverse_cum_sum = 0, direct_cum_sum = 0;
  for (int nn = larger_nn; nn <= larger_nn + smaller_nn; ++nn) {
    direct_cum_sum += distr[nn];
  }

  if (direct_cum_sum < tol) {
    direct_cum_sum = 1e10; // Any huge number
  }
  for (int nn = larger_nn + smaller_nn; nn > larger_nn; --nn) {
    double cur_val = distr[nn];
    reverse_cum_sum += cur_val;
    direct_cum_sum -= cur_val;

    res(nn - larger_nn, umi_ind) = reverse_cum_sum / direct_cum_sum;
  }
}

// [[Rcpp::export]]
NumericMatrix GetSmallerNeighboursDistributionsBySizes(const List &dp_matrices, const IntegerVector &larger_neighbours_num,
                                                       const s_vec_t &neighbour_prob_inds, int size_adj, int max_neighbour_num,
                                                       const IntegerVector &smaller_neighbours_num = IntegerVector()) {
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
    auto const &distr = matrix.column(size_adj);

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

    if (smaller_neighbours_num.size() != 0) {
      fillCumSumRatio(max_neighbour_num, smaller_neighbours_num.at(umi_ind), larger_nn, umi_ind, distr / prob_sum, res);
    }
    else {
      fillReverseCumSum(max_neighbour_num, larger_nn, umi_ind, distr / prob_sum, res);
    }
  }

  colnames(res) = as<StringVector>(larger_neighbours_num.names());

  return res;
}

// [[Rcpp::export]]
NumericVector GetSmallerNeighbourProbabilities(const NumericMatrix &small_neighs_dist, const IntegerVector &neighb_per_umi) {
  NumericVector res(sum(neighb_per_umi));
  int res_ind = 0;
  for (int umi_ind = 0; umi_ind < neighb_per_umi.size(); ++umi_ind) {
    int neighb_num = neighb_per_umi[umi_ind];
    for (int neighb_ind = 1; neighb_ind <= neighb_num; ++neighb_ind) {
      res[res_ind++] = small_neighs_dist(neighb_ind, umi_ind);
    }
  }

  return res;
}

// [[Rcpp::export]]
IntegerVector FilterUmisInGeneSimple(const IntegerVector &reads_per_umi, const std::vector<s_vec_t> &neighbourhood, double mult=1) {
  if (neighbourhood.size() != reads_per_umi.size())
    stop("Vectors must have equal size");

  if (reads_per_umi.size() == 1)
    return reads_per_umi;

  LogicalVector filt_mask(reads_per_umi.size(), true);

  Progress p(0, false);
  auto const &rpu_map = parseVector(reads_per_umi);
  auto const &umis = as<s_vec_t>(as<StringVector>(reads_per_umi.names()));

  std::vector<int> targets(umis.size(), -1);
  si_map_t umi_inds;
  for (int i = 0; i < umis.size(); ++i) {
    umi_inds[umis[i]] = i;
  }

  for (int i = 0; i < umis.size(); ++i) {
    auto const &neighbours = neighbourhood[i];
    if (neighbours.empty())
      continue;

    int cur_rpu = rpu_map.at(umis[i]);
    int neighb_index = -1;
    for (auto const &neighbour : neighbours) {
      auto neighb_it = rpu_map.find(neighbour);
      if (neighb_it == rpu_map.end() || (neighb_it->second < cur_rpu * mult - EPS)) // rpu_neighb >= cur_rpu
        continue;

      int cur_neighb_index = umi_inds.at(neighb_it->first);
      if (targets[cur_neighb_index] == i)
        continue;

      if (neighb_index != -1) {
        neighb_index = -1;
        break;
      }

      filt_mask[i] = false;
      neighb_index = cur_neighb_index;
    }

    if (neighb_index != -1) {
      targets[i] = neighb_index;
    }

    if (p.check_abort())
      break;
  }

  if (any(filt_mask).is_false()) {
    filt_mask[which_max(reads_per_umi)] = true;
  }

  return reads_per_umi[filt_mask];
}
