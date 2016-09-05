#include "PoissonRealBarcodesMergeStrategy.h"

#include <boost/iterator/transform_iterator.hpp>
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
		bool is_base_cb_real = base_cell_ind == neighbour_cells[0];
		double max_merge_prob = is_base_cb_real
								? PoissonRealBarcodesMergeStrategy::max_merge_prob
								: PoissonRealBarcodesMergeStrategy::max_real_cb_merge_prob;

		long best_target = -1;
		double min_prob = 2;
		for (auto cell_ind : neighbour_cells)
		{
			double prob = this->get_bootstrap_intersect_prob(container, base_cell_ind, cell_ind);

			container.stats().add(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
								  container.cell_barcode(cell_ind), 1);
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
																		  size_t fit_size, unsigned multiplies_count) const
	{
		auto const &cell1_dist = container.cell_genes(cell1_ind);
		auto const &cell2_dist = container.cell_genes(cell2_ind);
		size_t intersect_size = RealBarcodesMergeStrategy::get_umigs_intersect_size(cell1_dist, cell2_dist);
		if (intersect_size == 0)
			return 1;

		ul_list_t sizes;

		for (size_t repeats_count = 100; repeats_count < 1001 && repeats_count < fit_size; repeats_count += 900)
		{
			sizes.clear();
			sizes.reserve(repeats_count);
			double prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);

			if (prob > 0.1)
				return prob;
		}

		double estimated = -1, prob = -1;
		for (unsigned i = 1; i <= multiplies_count; ++i)
		{
			size_t repeats_count = fit_size * i; // TODO remove (it can fix only high probabilities, we don't care about it); replace 4 / size with std
			sizes.clear();
			sizes.reserve(repeats_count);
			prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);
			estimated = this->estimate_by_r(sizes, intersect_size);
			if (prob == 0 && estimated < 4. / repeats_count || estimated > 0 && prob > 0 && estimated / prob <= 4 && prob / estimated <= 4) // check for outliers
				return estimated;
		}

		L_WARN << "Not stable estimation: " << container.cell_barcode(cell1_ind) << ", " <<
					container.cell_barcode(cell2_ind) << "; probs: (" << estimated << ", " << prob << ")";
		return estimated;
	}

	double PoissonRealBarcodesMergeStrategy::estimate_by_r(ul_list_t sizes, size_t val) const
	{
		(*this->_r)["sizes"] = sizes;
		(*this->_r)["val"] = val - 1;
		this->_r->parseEvalQ("p_fit = fitdistr(sizes, \"poisson\")\n"
							  "res <- ppois(lambda=p_fit$estimate[1], q = val, lower.tail=F)");
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
			if (cell2_dist.find(item.first) == cell2_dist.end())
				continue;

			c1_dist.push_back(item.second.size());
			c2_dist.push_back(item.second.size());
		}

		size_t repeats_sum = 0;
		std::vector<unsigned> intersection_marks(this->_umis_number, 0);
		for (unsigned repeat_num = 1; repeat_num <= repeats_count; ++repeat_num)
		{
			size_t intersect_size = 0;
			auto cell2_it = c2_dist.begin();
			for (auto const &cell1_gene_size : c1_dist)
			{
				for (size_t choice_num = 0; choice_num < cell1_gene_size; ++choice_num)
				{
					intersection_marks[this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()]] = repeat_num;
				}

				size_t cell2_gene_size = *cell2_it;

				for (size_t choice_num = 0; choice_num < cell2_gene_size; ++choice_num)
				{
					bs_umi_t umi = this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()];
					if (intersection_marks[umi] == repeat_num)
					{
						intersect_size++;
						intersection_marks[umi] = 0;
					}
				}
				++cell2_it;
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