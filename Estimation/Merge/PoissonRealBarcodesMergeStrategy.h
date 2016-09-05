#pragma once

#include <RInside.h>
#include "RealBarcodesMergeStrategy.h"

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeProbs;
	struct testPoissonMergeInit;
	struct testPoissonMergeTime;
}

namespace Estimation
{
namespace Merge
{
	class PoissonRealBarcodesMergeStrategy : public RealBarcodesMergeStrategy
	{
		friend struct TestEstimatorMergeProbs::testPoissonMergeProbs;
		friend struct TestEstimatorMergeProbs::testPoissonMergeInit;
		friend struct TestEstimatorMergeProbs::testPoissonMergeTime;

	private:
		typedef Estimation::CellsDataContainer::s_ul_hash_t s_ul_hash_t;
		typedef unsigned bs_umi_t;

	private:
		static const double max_merge_prob;
		static const double max_real_cb_merge_prob;

		bs_umi_t _umis_number;

		s_ul_hash_t _umis_distribution;
		std::vector<bs_umi_t> _umis_bootstrap_distribution;
		RInside* _r;

	private:
		double get_bootstrap_intersect_prob(const CellsDataContainer &container, size_t cell1_ind, size_t cell2_ind,
											size_t fit_size = 3000, unsigned multiplies_count = 4) const;

	protected:
		virtual void init(const Estimation::CellsDataContainer &container) override;
		virtual void release() override;

		virtual long get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind, const ul_list_t &neighbour_cells) const override;
		virtual long get_max_merge_dist(long min_real_cb_dist) const override;

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