#pragma once

#include "MergeUMIsStrategyAbstract.h"

#include <unordered_set>
#include <unordered_map>

namespace TestEstimator
{
	struct testFillWrongUmis;
	struct testRemoveSimilarWrongUmis;
}

namespace Estimation
{
	namespace Merge
	{
		namespace UMIs
		{
			class MergeUMIsStrategySimple : public MergeUMIsStrategyAbstract
			{
				friend struct TestEstimator::testFillWrongUmis;
				friend struct TestEstimator::testRemoveSimilarWrongUmis;

			private:
				using s_hash_t = std::unordered_set<std::string>;
				using s_vec_t = std::vector<std::string>;

			private:
				const unsigned _max_merge_distance;
				static const std::string nucleotides;

			private:
				CellsDataContainer::s_s_hash_t find_targets(const StringIndexer &umi_indexer,
				                                            const Gene::umis_t &all_umis, const s_hash_t &bad_umis) const;
				void remove_similar_wrong_umis(s_vec_t &wrong_umis) const;
				CellsDataContainer::s_s_hash_t fill_wrong_umis(s_vec_t &wrong_umis) const;

				bool is_umi_real(const std::string &umi) const;

			public:
				explicit MergeUMIsStrategySimple(unsigned max_merge_distance);
				void merge(CellsDataContainer &container) const override;
			};
		}
	}
}
