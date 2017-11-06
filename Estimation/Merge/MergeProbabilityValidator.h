#pragma once

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
			typedef std::vector<size_t> ul_vec_t;
			typedef std::vector<double> d_vec_t;
		private:
			ul_vec_t _umis_per_cell1;
			ul_vec_t _umis_per_cell2;
			d_vec_t _merge_probs;
			ul_vec_t _intersection_size;
			d_vec_t _expected_intersection_size;

		public:
			void run_validation(const CellsDataContainer &container, unsigned min_ed, size_t cb_pairs_num);

			const ul_vec_t& umis_per_cell1() const;
			const ul_vec_t& umis_per_cell2() const;
			const d_vec_t& merge_probs() const;
			const ul_vec_t& intersection_size() const;
			const d_vec_t& expected_intersection_size() const;
		};
	}
}