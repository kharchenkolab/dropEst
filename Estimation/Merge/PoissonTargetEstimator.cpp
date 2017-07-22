#include <Tools/UtilFunctions.h>
#include "PoissonTargetEstimator.h"
#include "MergeStrategyBase.h"

namespace Estimation
{
namespace Merge
{
PoissonTargetEstimator::PoissonTargetEstimator(double max_merge_prob, double max_real_cb_merge_prob)
	: max_merge_prob(max_merge_prob)
	, max_real_cb_merge_prob(max_real_cb_merge_prob)
{}

long PoissonTargetEstimator::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
                                                   const ul_list_t &neighbour_cells) const
{
	bool is_base_cb_real = (base_cell_ind == neighbour_cells.at(0));
	double max_merge_prob = is_base_cb_real
							? PoissonTargetEstimator::max_merge_prob
							: PoissonTargetEstimator::max_real_cb_merge_prob;

	max_merge_prob /= neighbour_cells.size();

	long best_target = -1;
	double min_prob = 2;
	for (size_t cell_ind : neighbour_cells)
	{
		if (cell_ind == base_cell_ind)
			continue;

		double prob = this->get_bootstrap_intersect_prob(container, base_cell_ind, cell_ind);

//		container.stats().set(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
//							  container.cell_barcode(cell_ind), prob);
		if (prob < min_prob)
		{
			min_prob = prob;
			best_target = cell_ind;
		}
	}

	if (min_prob > max_merge_prob)
		return is_base_cb_real ? (long) base_cell_ind : -1;

	return best_target;
}

void PoissonTargetEstimator::init(const Estimation::CellsDataContainer &container)
{
	this->_umis_distribution = container.umis_distribution();
	bs_umi_t umi_ind = 0;
	for (auto const &umi: this->_umis_distribution)
	{
		this->_umis_bootstrap_distribution.insert(this->_umis_bootstrap_distribution.end(), umi.second, umi_ind);
		++umi_ind;
	}

	this->_umis_number = umi_ind;

	this->_r = Tools::init_r();
	this->_r->parseEvalQ("library(fitdistrplus)");

	srand(48);
}

void PoissonTargetEstimator::release()
{
	this->_umis_bootstrap_distribution.clear();
	this->_umis_distribution.clear();
}

double PoissonTargetEstimator::get_bootstrap_intersect_prob(const CellsDataContainer &container,
																	  size_t cell1_ind, size_t cell2_ind,
																	  size_t repeats_count, unsigned multiplies_count) const
{
	const size_t test_sample_size = 100;
	auto const &cell1_dist = container.cell_genes(cell1_ind);
	auto const &cell2_dist = container.cell_genes(cell2_ind);
	size_t intersect_size = MergeStrategyBase::get_umigs_intersect_size(cell1_dist, cell2_dist);
	if (intersect_size == 0)
		return 1;

	ul_list_t sizes;

	sizes.reserve(test_sample_size);
	double prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, test_sample_size, sizes);

	if (prob > 0.05) // Just to speed up
		return prob;

	double estimated = -1;
	for (unsigned i = 0; i < multiplies_count; ++i)
	{
		sizes.clear();
		sizes.reserve(repeats_count);
		this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);

		estimated = this->estimate_by_r(sizes, intersect_size);

		if (estimated >= 0)
			return estimated;

		repeats_count *= 2;
	}

	return -estimated;
}

double PoissonTargetEstimator::estimate_by_r(ul_list_t sizes, size_t val) const
{
	//TODO: optimize it with lambda = mean, "Use The d/q/p/q Statistical Functions (see RCpp book)"
	(*this->_r)["sizes"] = sizes;
	(*this->_r)["val"] = val;
	this->_r->parseEvalQ("p_fit <- fitdistr(sizes, \"poisson\")\n"
								 "res <- ppois(lambda=p_fit$estimate[1], q = val - 1, lower.tail=F)\n"
								 "res_upper <- max(ppois(lambda=p_fit$estimate[1] + 3 * p_fit$sd, q = val - 1, lower.tail=F),"
								 "ppois(lambda=max(p_fit$estimate[1] - 3 * p_fit$sd, 0), q = val - 1, lower.tail=F))\n"
								 "if (res_upper > 4 * res) {res <- -res_upper}\n"); // Check that we have enough observations
	return (*this->_r)["res"];
}

double PoissonTargetEstimator::get_bootstrap_intersect_sizes(const genes_t &cell1_dist,
																	   const genes_t &cell2_dist,
																	   size_t real_intersect_size,
																	   size_t repeats_count, ul_list_t &sizes) const
{
	std::vector<size_t >c1_dist, c2_dist;
	for (auto const &item : cell1_dist)
	{
		auto const& cell2_it = cell2_dist.find(item.first);
		if (cell2_it == cell2_dist.end())
			continue;

		c1_dist.push_back(item.second.size());
		c2_dist.push_back(cell2_it->second.size());
	}

	size_t repeats_sum = 0;
	std::vector<bs_umi_t> intersection_marks(this->_umis_number, 0);
	for (unsigned repeat_num = 1; repeat_num <= repeats_count; ++repeat_num)
	{
		size_t intersect_size = 0;
		for (auto c1_it = c1_dist.begin(), c2_it = c2_dist.begin(); c1_it != c1_dist.end(); ++c1_it, ++c2_it)
		{
			size_t cell1_gene_size = *c1_it;
			for (size_t choice_num = 0; choice_num < cell1_gene_size; ++choice_num)
			{
				intersection_marks[this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()]] = repeat_num;
			}

			size_t cell2_gene_size = *c2_it;
			for (size_t choice_num = 0; choice_num < cell2_gene_size; ++choice_num)
			{
				bs_umi_t umi = this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()];
				if (intersection_marks[umi] == repeat_num)
				{
					intersect_size++;
					intersection_marks[umi] = 0;
				}
			}
		}
		sizes.push_back(intersect_size);
		repeats_sum += intersect_size >= real_intersect_size;
	}

	return repeats_sum / (double) repeats_count;
}
}
}
