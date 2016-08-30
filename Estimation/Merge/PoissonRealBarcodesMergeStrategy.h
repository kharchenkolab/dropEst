#pragma once

#include <RInside.h>
#include "RealBarcodesMergeStrategy.h"

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeProbs;
	struct testPoissonMergeInit;
}

namespace Estimation
{
namespace Merge
{
	class PoissonRealBarcodesMergeStrategy : public RealBarcodesMergeStrategy
	{
		friend struct TestEstimatorMergeProbs::testPoissonMergeProbs;
		friend struct TestEstimatorMergeProbs::testPoissonMergeInit;

	private:
		typedef Estimation::CellsDataContainer::s_ul_hash_t s_ul_hash_t;

	private:
		static const double max_merge_prob;

		s_ul_hash_t _umis_distribution;
		names_t _umis_bootstrap_distribution;
		RInside* _r;

	private:
		double get_bootstrap_intersect_prob(const CellsDataContainer &container, size_t cell1_ind, size_t cell2_ind,
											size_t fit_size = 3000, unsigned multiplies_count = 4) const;

		static const s_ul_hash_t::key_type & get_key(const s_ul_hash_t::value_type & pair);
		static const names_t intersect_keys(const s_ul_hash_t &map1, const s_ul_hash_t &map2);

	protected:
		virtual long get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind, const ul_list_t &neighbour_cells) const override;

		virtual void init(const Estimation::CellsDataContainer &container) override;

		virtual void release() override;

	public:
		PoissonRealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
										 int min_genes_before_merge, int min_genes_after_merge,
										 int max_merge_edit_distance, double min_merge_fraction);

		double estimate_by_r(ul_list_t vector, size_t i) const;

		double get_bootstrap_intersect_sizes(const genes_t &cell1_dist, const genes_t &cell2_dist,
											 size_t real_intersect_size, size_t repeats_count, ul_list_t &sizes) const;
	};
}
}