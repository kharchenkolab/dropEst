#include "BrokenRealBarcodesMergeStrategy.h"

#include <Tools/UtilFunctions.h>
#include <Tools/Logs.h>

namespace Estimation
{
namespace Merge
{
	BrokenRealBarcodesMergeStrategy::BrokenRealBarcodesMergeStrategy(const std::string &barcodes_filename, size_t barcode2_length,
																	   int min_genes_before_merge, int min_genes_after_merge,
																	   int max_merge_edit_distance, double min_merge_fraction)
		: RealBarcodesMergeStrategy(barcodes_filename, barcode2_length, min_genes_before_merge,
									min_genes_after_merge, max_merge_edit_distance, min_merge_fraction)
	{ }

	long BrokenRealBarcodesMergeStrategy::get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind,
																 const ul_list_t &neighbour_cells) const
	{
		if (base_cell_ind == neighbour_cells[0])
			return base_cell_ind;

		return neighbour_cells[rand() % neighbour_cells.size()];
	}

	void BrokenRealBarcodesMergeStrategy::init(const Estimation::CellsDataContainer &container)
	{
		L_TRACE << "Broken merge selected";
		RealBarcodesMergeStrategy::init(container);
		srand(48);
	}

	long BrokenRealBarcodesMergeStrategy::get_max_merge_dist(long min_real_cb_dist) const
	{
		return min_real_cb_dist == 0 ? 0 : min_real_cb_dist + 15;
	}
}
}