#pragma once

#include <memory>
#include <vector>

#include "PoissonTargetEstimator.h"
#include "Estimation/CellsDataContainer.h"

namespace Estimation
{
	namespace Merge
	{
		class MergeProbabilityValidator
		{
		public:
			typedef std::vector<unsigned> u_vec_t;
			typedef std::vector<double> d_vec_t;

		private:
			u_vec_t _umis_per_cell1;
			u_vec_t _umis_per_cell2;
			u_vec_t _edit_distances;
			d_vec_t _merge_probs;
			u_vec_t _intersection_size;
			d_vec_t _expected_intersection_size;
			std::shared_ptr<PoissonTargetEstimator> _estimator;

		public:
			MergeProbabilityValidator(const std::shared_ptr<PoissonTargetEstimator> &_estimator);

			void run_validation(const CellsDataContainer &container, unsigned min_ed, unsigned max_ed,
			                    size_t cb_pairs_num, unsigned log_period);

			const u_vec_t& umis_per_cell1() const;
			const u_vec_t& umis_per_cell2() const;
			const u_vec_t& edit_distances() const;
			const d_vec_t& merge_probs() const;
			const u_vec_t& intersection_size() const;
			const d_vec_t& expected_intersection_size() const;
		};
	}
}