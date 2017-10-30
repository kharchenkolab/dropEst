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

long PoissonTargetEstimator::get_best_merge_target(CellsDataContainer &container, size_t base_cell_ind,
                                                   const ul_list_t &neighbour_cells)
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

		double prob = this->get_intersection_prob(container, base_cell_ind, cell_ind);

		container.cell(base_cell_ind).stats().set(Stats::MERGE_PROB_PER_TARGET_PER_CELL,
		                                          container.cell(cell_ind).barcode(), prob);

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

void PoissonTargetEstimator::init(const Estimation::CellsDataContainer::s_ul_hash_t &umi_distribution)
{
	double sum = 0;
	for (auto const &it : umi_distribution)
	{
		sum += it.second;
	}

	for (auto const &it : umi_distribution)
	{
		this->_umi_distribution.push_back(it.second / sum);
	}

	this->_umi_sampler = std::discrete_distribution<umi_t>(this->_umi_distribution.begin(), this->_umi_distribution.end());
//	this->_adjuster.init(this->_umi_distribution);

	Tools::init_r();
}

void PoissonTargetEstimator::release()
{
	this->_umi_distribution.clear();
}

double PoissonTargetEstimator::get_intersection_prob(const CellsDataContainer &container, size_t cell1_ind,
                                                     size_t cell2_ind)
{
	auto const &cell1 = container.cell(cell1_ind);
	auto const &cell2 = container.cell(cell2_ind);
	size_t intersect_size = MergeStrategyBase::get_umigs_intersect_size(cell1, cell2);
	if (intersect_size == 0)
		return 1;

	double est_intersection_size = 0;
	for (auto const &gene1_it : cell1.genes())
	{
		auto const& gene2_it = cell2.genes().find(gene1_it.first);
		if (gene2_it == cell2.genes().end())
			continue;

		est_intersection_size += this->estimate_genes_intersection_size(gene1_it.second.size(), gene2_it->second.size());
	}

	return Rcpp::ppois(Rcpp::IntegerVector::create(intersect_size - 1), est_intersection_size, false)[0];
}

double PoissonTargetEstimator::estimate_genes_intersection_size(size_t gene1_size, size_t gene2_size)
{
	const size_t repeats_num = 10000;

	if (gene1_size > gene2_size)
	{
		std::swap(gene1_size, gene2_size);
	}

	auto gene_pair = std::make_pair(gene1_size, gene2_size);
	auto pair_it = this->_estimated_gene_intersections.find(gene_pair);
	if (pair_it != this->_estimated_gene_intersections.end())
		return pair_it->second;

	size_t repeats_sum = 0;
	std::vector<long> intersection_marks(this->_umi_distribution.size(), 0);
	std::mt19937 gen;

	for (unsigned repeat_id = 1; repeat_id <= repeats_num; ++repeat_id)
	{
		size_t intersect_size = 0;

		for (size_t choice_num = 0; choice_num < gene1_size; ++choice_num)
		{
			umi_t sample_umi = this->_umi_sampler(gen);
			if (intersection_marks[sample_umi] == repeat_id)
			{
				--choice_num;
				continue;
			}

			intersection_marks[sample_umi] = repeat_id;
		}

		for (size_t choice_num = 0; choice_num < gene2_size; ++choice_num)
		{
			umi_t sample_umi = this->_umi_sampler(gen);
			if (intersection_marks[sample_umi] == -repeat_id)
			{
				--choice_num;
				continue;
			}

			if (intersection_marks[sample_umi] == repeat_id)
			{
				intersect_size++;
			}

			intersection_marks[sample_umi] = -repeat_id;
		}

		repeats_sum += intersect_size;
	}

	double est_size = repeats_sum / (double) repeats_num;

	this->_estimated_gene_intersections.emplace(gene_pair, est_size);
	return est_size;
}
}
}
