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
			virtual long get_merge_target(CellsDataContainer &container, size_t base_cell_ind) override;
			virtual size_t get_log_period() const override;

		public:
			PoissonSimpleMergeStrategy(const PoissonTargetEstimator &target_estimator, unsigned min_genes_before_merge,
			                           unsigned min_genes_after_merge, unsigned max_merge_edit_distance);

			virtual std::string merge_type() const override;
		};
	}
}