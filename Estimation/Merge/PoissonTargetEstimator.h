#pragma once


#include <Estimation/Cell.h>
#include <Estimation/CellsDataContainer.h>
#include <RInside.h>

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeProbs;
	struct testPoissonMergeInit;
	struct testPoissonMergeTime;
	struct testPoissonMergeRejections;
}

namespace Estimation
{
	namespace Merge
	{
		class PoissonTargetEstimator
		{
			friend struct TestEstimatorMergeProbs::testPoissonMergeProbs;
			friend struct TestEstimatorMergeProbs::testPoissonMergeInit;
			friend struct TestEstimatorMergeProbs::testPoissonMergeTime;

		private:
			typedef Estimation::CellsDataContainer::s_ul_hash_t s_ul_hash_t;
			typedef Estimation::CellsDataContainer::counts_t ul_list_t;
			typedef unsigned bs_umi_t;

		private:
			const double max_merge_prob;
			const double max_real_cb_merge_prob;

			bs_umi_t _umis_number;

			s_ul_hash_t _umis_distribution;
			std::vector<bs_umi_t> _umis_bootstrap_distribution;
			RInside* _r;

		private:
			double get_bootstrap_intersect_prob(const CellsDataContainer &container, size_t cell1_ind, size_t cell2_ind,
												size_t repeats_count = 700, unsigned multiplies_count = 3) const;

			double estimate_by_r(ul_list_t vector, size_t i) const;

			double get_bootstrap_intersect_sizes(const Cell &cell1, const Cell &cell2, size_t real_intersect_size,
			                                     size_t repeats_count, ul_list_t &sizes) const;

		public:
			PoissonTargetEstimator(double max_merge_prob, double max_real_cb_merge_prob);
			virtual void init(const Estimation::CellsDataContainer &container);
			virtual void release();

			virtual long get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
											   const ul_list_t &neighbour_cells) const;
		};
	}
}

