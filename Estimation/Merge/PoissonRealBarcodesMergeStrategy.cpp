#include "PoissonRealBarcodesMergeStrategy.h"

#include <Tools/UtilFunctions.h>
#include <Tools/Logs.h>

namespace Estimation
{
namespace Merge
{
	const double PoissonRealBarcodesMergeStrategy::max_merge_prob = 1e-5;
	const double PoissonRealBarcodesMergeStrategy::max_real_cb_merge_prob = 1e-7;

	PoissonRealBarcodesMergeStrategy::PoissonRealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
																	   int min_genes_before_merge, int min_genes_after_merge,
																	   int max_merge_edit_distance, double min_merge_fraction)
		: RealBarcodesMergeStrategy(barcodes_filename, barcode2_length, min_genes_before_merge,
									min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
	{ }

	long PoissonRealBarcodesMergeStrategy::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
																 const ul_list_t &neighbour_cells) const
	{
		bool is_base_cb_real = (base_cell_ind == neighbour_cells[0]);
		double max_merge_prob = is_base_cb_real
								? PoissonRealBarcodesMergeStrategy::max_merge_prob
								: PoissonRealBarcodesMergeStrategy::max_real_cb_merge_prob;

		max_merge_prob /= neighbour_cells.size();

		long best_target = -1;
		double min_prob = 2;
		for (auto cell_ind : neighbour_cells)
		{
			if (cell_ind == base_cell_ind)
				continue;

			double prob = this->get_bootstrap_intersect_prob(container, base_cell_ind, cell_ind);

			container.stats().set(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
								  container.cell_barcode(cell_ind), prob);
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

	void PoissonRealBarcodesMergeStrategy::init(const Estimation::CellsDataContainer &container)
	{
		RealBarcodesMergeStrategy::init(container);
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

	void PoissonRealBarcodesMergeStrategy::release()
	{
		RealBarcodesMergeStrategy::release();
		this->_umis_bootstrap_distribution.clear();
		this->_umis_distribution.clear();
	}

	double PoissonRealBarcodesMergeStrategy::get_bootstrap_intersect_prob(const CellsDataContainer &container,
																		  size_t cell1_ind, size_t cell2_ind,
																		  size_t repeats_count, unsigned multiplies_count) const
	{
		const size_t test_sample_size = 100;
		auto const &cell1_dist = container.cell_genes(cell1_ind);
		auto const &cell2_dist = container.cell_genes(cell2_ind);
		size_t intersect_size = RealBarcodesMergeStrategy::get_umigs_intersect_size(cell1_dist, cell2_dist);
		if (intersect_size == 0)
			return 1;

		ul_list_t sizes;

		sizes.reserve(test_sample_size);
		double prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, test_sample_size, sizes);

		if (prob > 0.05) // Just to speed up
			return prob;

		for (unsigned i = 0; i < multiplies_count; ++i)
		{
			sizes.clear();
			sizes.reserve(repeats_count);
			this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);
			double estimated = this->estimate_by_r(sizes, intersect_size);

			if (estimated != -1)
				return estimated;

			repeats_count *= 2;
		}

		return 2;
	}

	double PoissonRealBarcodesMergeStrategy::estimate_by_r(ul_list_t sizes, size_t val) const
	{
		(*this->_r)["sizes"] = sizes;
		(*this->_r)["val"] = val;
		this->_r->parseEvalQ("p_fit <- fitdistr(sizes, \"poisson\")\n"
							  "res <- ppois(lambda=p_fit$estimate[1], q = val - 1, lower.tail=F)\n"
							  "res_upper <- ppois(lambda=p_fit$estimate[1] + 3 * p_fit$sd, q = val - 1, lower.tail=F)\n"
							  "if (res_upper > 4 * res) {res <- -1}\n"); // Check that we have enough observations
		return (*this->_r)["res"];
	}

	double PoissonRealBarcodesMergeStrategy::get_bootstrap_intersect_sizes(const genes_t &cell1_dist,
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

	long PoissonRealBarcodesMergeStrategy::get_max_merge_dist(long min_real_cb_dist) const
	{
		return min_real_cb_dist == 0 ? 2 : min_real_cb_dist + 1;
	}
}
}