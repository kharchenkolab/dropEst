#pragma once

#include <Estimation/CellsDataContainer.h>
#include "AbstractMergeStrategy.h"

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
	class RealBarcodesMergeStrategy : public AbstractMergeStrategy
	{
		friend struct TestEstimator::testBarcodesFile;
		friend struct TestEstimator::testUmigsIntersection;
		friend struct TestEstimator::testFillDistances;
		friend struct TestEstimator::testRealNeighboursCbs;
		friend struct TestEstimator::testRealNeighbours;

	private:
		std::string _barcodes_filename;
		size_t _barcode2_length;

		static const int MAX_REAL_MERGE_EDIT_DISTANCE=5;

	private:
		ids_t get_real_neighbour_cbs(const Estimation::CellsDataContainer &container, const names_t &cbs1, const names_t &cbs2,
									 const std::string &base_cb, const i_counter_t &dists1, const i_counter_t &dists2) const;

		size_t get_real_cb(const Estimation::CellsDataContainer &container, size_t base_cell_ind, const names_t &cbs1,
						 const names_t &cbs2) const;

		double get_umigs_intersect_fraction(const Estimation::CellsDataContainer &container, size_t cell1_ind,
											size_t cell2_ind) const;

		static void get_barcodes_list(const std::string &barcodes_filename, std::vector<std::string> &barcodes1,
									  std::vector<std::string> &barcodes2);

		static void fill_distances_to_cb(const std::string &cb_part1, const std::string &cb_part2, const names_t &cbs1,
										 const names_t &cbs2, i_counter_t &dists1, i_counter_t &dists2);

	public:
		RealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
								  int min_genes_before_merge, int min_genes_after_merge,
								  int max_merge_edit_distance, double min_merge_fraction);

		virtual void merge(Estimation::CellsDataContainer &container, const s_ii_hash_t &umig_cells_counts,
						   ids_t &filtered_cells) const override;
	};
}


}