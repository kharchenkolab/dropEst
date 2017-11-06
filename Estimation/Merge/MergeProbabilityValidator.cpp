#include <Tools/Logs.h>
#include "MergeProbabilityValidator.h"

namespace Estimation
{
	namespace Merge
	{
		void MergeProbabilityValidator::run_validation(const CellsDataContainer &container, unsigned min_ed,
		                                               unsigned max_ed, size_t cb_pairs_num)
		{
			if (container.filtered_cells().empty())
				return;

			srand(42);
			const size_t max_cell_id = container.filtered_cells().size();
			PoissonTargetEstimator estimator(1, 1);
			estimator.init(container.umi_distribution());
			for (size_t iter_num = 0; iter_num < cb_pairs_num; ++iter_num)
			{
				if (iter_num % 100000 == 0)
				{
					Tools::trace_time("Iteration: " + std::to_string(iter_num) + ": cache size " + std::to_string(estimator.cache_size()));
				}

				size_t cb_id1, cb_id2;
				while (true)
				{
					cb_id1 = container.filtered_cells()[rand() % max_cell_id];
					cb_id2  = container.filtered_cells()[rand() % max_cell_id];
					if (cb_id1 == cb_id2)
						continue;

					auto edit_distance = Tools::edit_distance(container.cell(cb_id1).barcode_c(), container.cell(cb_id2).barcode_c(), true, min_ed);
					if (edit_distance >= min_ed && edit_distance <= max_ed)
						break;
				}

				auto estimation_result = estimator.estimate_intersection_prob(container, cb_id1, cb_id2);

				this->_umis_per_cell1.push_back(container.cell(cb_id1).umis_number());
				this->_umis_per_cell2.push_back(container.cell(cb_id2).umis_number());
				this->_merge_probs.push_back(estimation_result.merge_probability);
				this->_intersection_size.push_back(estimation_result.intersection_size);
				this->_expected_intersection_size.push_back(estimation_result.expected_intersection_size);
			}
		}

		const MergeProbabilityValidator::ul_vec_t &MergeProbabilityValidator::umis_per_cell1() const
		{
			return this->_umis_per_cell1;
		}

		const MergeProbabilityValidator::ul_vec_t &MergeProbabilityValidator::umis_per_cell2() const
		{
			return this->_umis_per_cell2;
		}

		const MergeProbabilityValidator::d_vec_t &MergeProbabilityValidator::merge_probs() const
		{
			return this->_merge_probs;
		}

		const MergeProbabilityValidator::ul_vec_t &MergeProbabilityValidator::intersection_size() const
		{
			return this->_intersection_size;
		}

		const MergeProbabilityValidator::d_vec_t &MergeProbabilityValidator::expected_intersection_size() const
		{
			return this->_expected_intersection_size;
		}
	}
}