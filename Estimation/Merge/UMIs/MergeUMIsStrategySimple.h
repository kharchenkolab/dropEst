#pragma once

#include <Estimation/Cell.h>
#include <Estimation/CellsDataContainer.h>

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
			class MergeUMIsStrategySimple
			{
				friend struct TestEstimator::testFillWrongUmis;
				friend struct TestEstimator::testRemoveSimilarWrongUmis;

			private:
				typedef std::unordered_set<std::string> s_hash_t;
				typedef std::vector<std::string> s_vec_t;
			private:
				const unsigned _max_merge_distance;
				static const std::string nucleotides;

			private:
				CellsDataContainer::s_s_hash_t find_targets(const Cell::umi_map_t &all_umis, const s_hash_t &bad_umis) const;
				void remove_similar_wrong_umis(s_vec_t &wrong_umis) const;
				CellsDataContainer::s_s_hash_t fill_wrong_umis(s_vec_t &wrong_umis) const;

			public:
				MergeUMIsStrategySimple(unsigned max_merge_distance);
				void merge(CellsDataContainer &container) const;

				bool is_umi_real(const std::string &umi) const;
			};
		}
	}
}