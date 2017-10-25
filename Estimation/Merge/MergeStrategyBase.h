#pragma once

#include "MergeStrategyAbstract.h"
#include <Estimation/Cell.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <cstdlib>

namespace TestEstimator
{
	struct testMerge;
}

namespace Estimation
{
namespace Merge
{
	class MergeStrategyBase : public MergeStrategyAbstract
	{
		friend struct TestEstimator::testMerge;

	protected:
		typedef boost::unordered_map<size_t, size_t> u_u_hash_t;
		typedef boost::unordered_set<size_t> id_set_t;
		typedef boost::unordered_map<size_t, id_set_t> id_id_set_map_t;

	protected:
		const double _min_merge_fraction;
		const int _max_merge_edit_distance;

	private:
		void reassign(size_t cell_id, size_t target_cell_id, ul_list_t &cb_reassign_targets,
					  id_id_set_map_t &cb_reassigned_to_it) const;

	protected:
		void merge_force(Estimation::CellsDataContainer &container, size_t src_cell_id, size_t target_cell_ind,
						 ul_list_t &cb_reassign_targets, id_id_set_map_t &cb_reassigned_to_it) const;

		virtual ul_list_t merge_inited(Estimation::CellsDataContainer &container) override;
		virtual long get_merge_target(CellsDataContainer &container, size_t base_cell_ind) = 0;

		virtual size_t get_log_period() const;

	public:
		MergeStrategyBase(size_t min_genes_before_merge, size_t min_genes_after_merge,
		                  unsigned max_merge_edit_distance, double min_merge_fraction);

		static size_t get_umigs_intersect_size(const Cell &cell1, const Cell &cell2);
	};
}
}