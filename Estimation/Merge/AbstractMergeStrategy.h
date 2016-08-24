#pragma once

#include "IMergeStrategy.h"

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
	class AbstractMergeStrategy : public IMergeStrategy
	{
		friend struct TestEstimator::testMerge;

	protected:
		typedef boost::unordered_set<size_t> i_set_t;
		typedef boost::unordered_map<size_t, i_set_t> ISIHM;

	protected:
		const double _min_merge_fraction;
		const int _max_merge_edit_distance;

	public:
		AbstractMergeStrategy(int min_genes_before_merge, int min_genes_after_merge, int max_merge_edit_distance,
							  double min_merge_fraction);

	private:
		void reassign(size_t cell_id, size_t target_cell_id, ids_t &cb_reassign_targets,
					  ISIHM &cb_reassigned_to_it) const;

	protected:
		void merge_force(Estimation::CellsDataContainer &container, size_t src_cell_id, size_t target_cell_ind,
						 int target_genes_count, ids_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it) const;
	};
}
}