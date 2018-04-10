#pragma once

#include "MergeUMIsStrategyAbstract.h"

namespace TestEstimator
{
	struct testUMIMergeStrategyDirectional;
}

namespace Estimation
{
	namespace Merge
	{
		namespace UMIs
		{
			class MergeUMIsStrategyDirectional : public MergeUMIsStrategyAbstract
			{
				friend struct TestEstimator::testUMIMergeStrategyDirectional;

			private:
				struct UmiWrap
				{
					std::string sequence;
					size_t n_reads;

					UmiWrap(const std::string &sequence, size_t n_reads);
				};

				using umi_vec_t = std::vector<UmiWrap>;
				using merge_targets_t = CellsDataContainer::s_s_hash_t;

			private:
				const double _mult;
				const unsigned _max_edit_distance;

			private:
				merge_targets_t find_targets(umi_vec_t &umis) const;

			public:
				explicit MergeUMIsStrategyDirectional(double mult = 2, unsigned max_edit_distance=1);
				void merge(CellsDataContainer &container) const override;
			};
		}
	}
}