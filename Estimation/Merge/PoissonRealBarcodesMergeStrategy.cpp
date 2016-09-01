#include "PoissonRealBarcodesMergeStrategy.h"

#include <boost/iterator/transform_iterator.hpp>
#include <Tools/UtilFunctions.h>
#include <Tools/Logs.h>

namespace Estimation
{
namespace Merge
{
	const double PoissonRealBarcodesMergeStrategy::max_merge_prob = 1e-5;

	const PoissonRealBarcodesMergeStrategy::s_ul_hash_t::key_type & PoissonRealBarcodesMergeStrategy::get_key(const s_ul_hash_t::value_type & pair)
	{
		return pair.first;
	}

	PoissonRealBarcodesMergeStrategy::PoissonRealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
																	   int min_genes_before_merge, int min_genes_after_merge,
																	   int max_merge_edit_distance, double min_merge_fraction)
		: RealBarcodesMergeStrategy(barcodes_filename, barcode2_length, min_genes_before_merge,
									min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
	{ }

	long PoissonRealBarcodesMergeStrategy::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
																 const ul_list_t &neighbour_cells) const
	{
		long best_target = -1;
		double min_prob = 2;
		for (auto cell_ind : neighbour_cells)
		{
			double prob = this->get_bootstrap_intersect_prob(container, base_cell_ind, cell_ind);

			container.stats().add(Stats::MERGE_PROB_BY_CELL, container.cell_barcode(base_cell_ind),
								  container.cell_barcode(cell_ind), 1);
			if (prob > min_prob)
			{
				min_prob = prob;
				best_target = cell_ind;
			}
			//TODO Tests
			//TODO Replace strings by ints in the bootstrap array
		}

		if (min_prob > PoissonRealBarcodesMergeStrategy::max_merge_prob)
			return -1;

		return best_target;
	}

	void PoissonRealBarcodesMergeStrategy::init(const Estimation::CellsDataContainer &container)
	{
		RealBarcodesMergeStrategy::init(container);
		this->_umis_distribution = container.umis_distribution();
		int umi_ind = 0;
		for (auto const &umi: this->_umis_distribution)
		{
			this->_umis_bootstrap_distribution.insert(this->_umis_bootstrap_distribution.end(), umi.second, umi_ind);
			++umi_ind;
		}

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
			size_t repeats_count = fit_size * i; // TODO remove (it can fix only high probabilities, we don't care about it)
			sizes.clear();
			sizes.reserve(repeats_count);
			prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);
			estimated = this->estimate_by_r(sizes, intersect_size);
			if (prob == 0 && estimated < 2 / sizes.size() || estimated > 0 && prob > 0 && estimated / prob <= 4 && prob / estimated <= 4) // check for outliers
				return estimated;
		}

		L_WARN << "Not stable estimation: " << container.cell_barcode(cell1_ind) << ", " <<
					container.cell_barcode(cell2_ind) << "; probs: (" << estimated << ", " << prob << ")";
		return estimated;
	}

	const PoissonRealBarcodesMergeStrategy::names_t PoissonRealBarcodesMergeStrategy::intersect_keys(const s_ul_hash_t &map1, const s_ul_hash_t &map2)
	{
		typedef const s_ul_hash_t::key_type & (*get_key_t)(const s_ul_hash_t::value_type &);
		typedef boost::transform_iterator<get_key_t, s_ul_hash_t::const_iterator> key_iterator_t;

		names_t intersection;
		std::set_intersection(key_iterator_t(map1.begin(), PoissonRealBarcodesMergeStrategy::get_key),
							  key_iterator_t(map1.end(), PoissonRealBarcodesMergeStrategy::get_key),
							  key_iterator_t(map2.begin(), PoissonRealBarcodesMergeStrategy::get_key),
							  key_iterator_t(map2.end(), PoissonRealBarcodesMergeStrategy::get_key),
							  std::back_inserter(intersection));

		return intersection;
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
		size_t repeats_sum = 0;
		for (int repeat_num = 0; repeat_num < repeats_count; ++repeat_num) //TODO Optimize (common genes)
		{
			size_t intersect_size = 0;
			for (auto const &gene : cell1_dist)
			{
				std::set<bs_umi_t> umis_set;
				size_t cell1_gene_size = cell1_dist.at(gene.first).size();
				auto cell2_it = cell2_dist.find(gene.first);
				if (cell2_it == cell2_dist.end())
					continue;

				for (size_t choice_num = 0; choice_num < cell1_gene_size; ++choice_num)
				{
					umis_set.emplace(this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()]);
				}

				size_t cell2_gene_size = cell2_it->second.size();
				for (size_t choice_num = 0; choice_num < cell2_gene_size; ++choice_num)
				{
					const bs_umi_t &umi = this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()];
					if (umis_set.find(umi) != umis_set.end())
					{
						intersect_size++;
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