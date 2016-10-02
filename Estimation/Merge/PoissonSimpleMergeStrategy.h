#pragma once


#include "SimpleMergeStrategy.h"
#include "PoissonTargetEstimator.h"

namespace Estimation
{
	namespace Merge
	{
		class PoissonSimpleMergeStrategy : public SimpleMergeStrategy
		{
		private:
			PoissonTargetEstimator _target_estimator;

		protected:
			virtual void init(const CellsDataContainer &container) override;
			virtual void release() override;

		public:
			PoissonSimpleMergeStrategy(int min_genes_before_merge, int min_genes_after_merge,
									   int max_merge_edit_distance, double min_merge_fraction);

			virtual long get_merge_target(const CellsDataContainer &container, size_t base_cell_ind) const override;

			virtual std::string merge_type() const override;
		};
	}
}