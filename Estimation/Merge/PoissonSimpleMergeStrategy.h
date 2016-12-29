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
			PoissonSimpleMergeStrategy(const boost::property_tree::ptree &config);

			virtual long get_merge_target(const CellsDataContainer &container, size_t base_cell_ind) const override;

			virtual std::string merge_type() const override;
		};
	}
}