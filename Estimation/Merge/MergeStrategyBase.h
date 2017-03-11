#pragma once

#include "MergeStrategyAbstract.h"

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
		typedef boost::unordered_set<size_t> i_set_t;
		typedef boost::unordered_map<size_t, i_set_t> ISIHM;
		typedef Estimation::CellsDataContainer::genes_t genes_t;

	protected:
		const double _min_merge_fraction;
		const int _max_merge_edit_distance;

	private:
		void reassign(size_t cell_id, size_t target_cell_id, ul_list_t &cb_reassign_targets,
					  ISIHM &cb_reassigned_to_it) const;

	protected:
		void merge_force(Estimation::CellsDataContainer &container, size_t src_cell_id, size_t target_cell_ind,
						 ul_list_t &cb_reassign_targets, ISIHM &cb_reassigned_to_it) const;

		virtual ul_list_t merge_inited(Estimation::CellsDataContainer &container) const override;
		virtual long get_merge_target(const Estimation::CellsDataContainer &container, size_t base_cell_ind) const = 0;

		virtual size_t get_log_period() const;

	public:
		MergeStrategyBase(const boost::property_tree::ptree &config);

		static size_t get_umigs_intersect_size(const genes_t &cell1_dist, const genes_t &cell2_dist);
	};
}
}