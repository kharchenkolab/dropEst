#include "PoissonRealBarcodesMergeStrategy.h"

#include <boost/iterator/transform_iterator.hpp>
#include <RInside.h>
#include <Tools/UtilFunctions.h>

namespace Estimation
{
namespace Merge
{
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
		double min_prob = 1;
		for (auto cell_ind : neighbour_cells)
		{
			//TODO call get_bootstrap_intersect_prob
			//TODO check if prob == -1
			//TODO check if min prob is to high
			//TODO Tests
		}

		return RealBarcodesMergeStrategy::get_best_merge_target(container, base_cell_ind, neighbour_cells);
	}

	void PoissonRealBarcodesMergeStrategy::init(const Estimation::CellsDataContainer &container)
	{
		RealBarcodesMergeStrategy::init(container);
		this->_umis_distribution = container.umis_distribution();
		for (auto const &umi: this->_umis_distribution)
		{
			this->_umis_bootstrap_distribution.insert(this->_umis_bootstrap_distribution.end(), umi.second, umi.first);
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

	double PoissonRealBarcodesMergeStrategy::get_bootstrap_intersect_prob(const genes_t &cell1_dist,
																		  const genes_t &cell2_dist,
																		  size_t intersect_size, size_t fit_size,
																		  unsigned multiplies_count)
	{
		if (intersect_size == 0)
			return 1;

		ul_list_t sizes;

		for (size_t repeats_count = 100; repeats_count < 1001 && repeats_count < fit_size; repeats_count += 900)
		{
			sizes.clear();
			sizes.reserve(repeats_count);
			double prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);

			if (prob > 0.05)
				return prob;
		}

		double estimated = -1;
		for (unsigned i = 1; i <= multiplies_count; ++i)
		{
			size_t repeats_count = fit_size * i;
			sizes.clear();
			sizes.reserve(repeats_count);
			double prob = this->get_bootstrap_intersect_sizes(cell1_dist, cell2_dist, intersect_size, repeats_count, sizes);
			estimated = this->estimate_by_r(sizes, intersect_size);
			if (prob == 0 && estimated < 2 / sizes.size() || estimated / prob <= 4 && prob / estimated <= 4) // check for outliers
				return estimated;
		}

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

	double PoissonRealBarcodesMergeStrategy::estimate_by_r(ul_list_t sizes, size_t val)
	{
		(*this->_r)["sizes"] = sizes;
		(*this->_r)["val"] = val;
		this->_r->parseEvalQ("p_fit = fitdistr(sizes, \"poisson\")\n"
							  "res <- dpois(lambda=p_fit$estimate[1], x = val)");
		return (*this->_r)["res"];
	}

	double PoissonRealBarcodesMergeStrategy::get_bootstrap_intersect_sizes(const genes_t &cell1_dist,
																		   const genes_t &cell2_dist,
																		   size_t real_intersect_size,
																		   size_t repeats_count, ul_list_t sizes)
	{
//		names_t common_genes(PoissonRealBarcodesMergeStrategy::intersect_keys(cell1_dist, cell2_dist));
		size_t repeats_sum = 0;
		for (int repeat_num = 0; repeat_num < repeats_count; ++repeat_num)
		{
			std::set<std::string> gene_set;
			size_t intersect_size = 0;
			for (auto const &gene : cell1_dist)
			{
				size_t gene1_size = cell1_dist.at(gene.first).size();
				auto cell2_it = cell2_dist.find(gene.first);
				if (cell2_it == cell2_dist.end())
					continue;

				for (int choice_num = 0; choice_num < gene1_size; ++choice_num)
				{
					gene_set.emplace(this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()]);
				}

				size_t gene2_size = cell2_it->second.size();
				for (int choice_num = 0; choice_num < gene2_size; ++choice_num)
				{
					const std::string &umi = this->_umis_bootstrap_distribution[rand() % this->_umis_bootstrap_distribution.size()];
					if (gene_set.find(umi) != gene_set.end())
					{
						intersect_size++;
					}
				}
			}
			sizes.push_back(intersect_size);
			repeats_sum += intersect_size >= real_intersect_size;
		}

		return repeats_sum / repeats_count;
	}
}
}