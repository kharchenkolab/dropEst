#pragma once

#include <Estimation/CellsDataContainer.h>
#include "MergeStrategyBase.h"

namespace TestEstimator
{
	struct testBarcodesFile;
	struct testUmigsIntersection;
	struct testFillDistances;
	struct testRealNeighboursCbs;
	struct testRealNeighbours;
}

namespace Estimation
{
namespace Merge
{
	class RealBarcodesMergeStrategy : public MergeStrategyBase
	{
		friend struct TestEstimator::testBarcodesFile;
		friend struct TestEstimator::testUmigsIntersection;
		friend struct TestEstimator::testFillDistances;
		friend struct TestEstimator::testRealNeighboursCbs;
		friend struct TestEstimator::testRealNeighbours;

	protected:
		typedef Estimation::CellsDataContainer::genes_t genes_t;

	private:
		std::string _barcodes_filename;
		size_t _barcode2_length;

		names_t _barcodes1;
		names_t _barcodes2;

		static const int MAX_REAL_MERGE_EDIT_DISTANCE=5;

	private:
		ul_list_t get_real_neighbour_cbs(const Estimation::CellsDataContainer &container, const std::string &base_cb,
									 const i_counter_t &dists1, const i_counter_t &dists2) const;

		long get_real_cb(const Estimation::CellsDataContainer &container, size_t base_cell_ind) const;

		static double get_umigs_intersect_fraction(const genes_t &cell1_dist, const genes_t &cell2_dist);

		static void get_barcodes_list(const std::string &barcodes_filename, std::vector<std::string> &barcodes1,
									  std::vector<std::string> &barcodes2);

		void fill_distances_to_cb(const std::string &cb_part1, const std::string &cb_part2, i_counter_t &dists1,
								  i_counter_t &dists2) const;

	protected:
		virtual void init(const Estimation::CellsDataContainer &container) override;
		virtual void release() override;

		virtual long get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind, const ul_list_t &neighbour_cells) const;
		virtual long get_max_merge_dist(long min_real_cb_dist) const;

		virtual bool need_parallel() const;

		static size_t get_umigs_intersect_size(const genes_t &cell1_dist, const genes_t &cell2_dist);

	public:
		RealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
								  int min_genes_before_merge, int min_genes_after_merge,
								  int max_merge_edit_distance, double min_merge_fraction);

		virtual void merge_inited(Estimation::CellsDataContainer &container, const s_uu_hash_t &umig_cells_counts,
						   ul_list_t &filtered_cells) const override;
	};
}


}