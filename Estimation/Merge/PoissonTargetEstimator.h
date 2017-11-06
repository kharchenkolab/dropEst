#pragma once


#include <Estimation/Cell.h>
#include <Estimation/CellsDataContainer.h>

#include <Tools/UtilFunctions.h>
#include <Tools/CollisionsAdjuster.h>
#include <random>

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeProbs;
	struct testPoissonMergeInit;
	struct testPoissonMergeTime;
	struct testPoissonMergeRejections;
	struct testIntersectionSizeEstimation;
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
			friend struct TestEstimatorMergeProbs::testIntersectionSizeEstimation;

		public:
			class EstimationResult
			{
			public:
				const size_t intersection_size;
				const double expected_intersection_size;
				const double merge_probability;

			public:
				EstimationResult(size_t intersection_size, double expected_intersection_size, double merge_probability);
			};

		private:
			typedef Estimation::CellsDataContainer::ids_t ul_list_t;
			typedef unsigned umi_t;

		private:
			const double max_merge_prob;
			const double max_real_cb_merge_prob;

			std::discrete_distribution<umi_t> _umi_sampler;
//			Tools::CollisionsAdjuster _adjuster;
			std::vector<double> _umi_distribution;
			std::unordered_map<std::pair<size_t, size_t>, double, Tools::PairHash> _estimated_gene_intersections;

		private:
			double estimate_genes_intersection_size(size_t gene1_size, size_t gene2_size);

		public:
			PoissonTargetEstimator(double max_merge_prob, double max_real_cb_merge_prob);
			virtual void init(const Estimation::CellsDataContainer::s_ul_hash_t &umi_distribution);
			virtual void release();

			size_t cache_size() const;
			EstimationResult estimate_intersection_prob(const CellsDataContainer &container, size_t cell1_ind, size_t cell2_ind);
			virtual long get_best_merge_target(CellsDataContainer &container, size_t base_cell_ind,
			                                   const ul_list_t &neighbour_cells);
		};
	}
}

