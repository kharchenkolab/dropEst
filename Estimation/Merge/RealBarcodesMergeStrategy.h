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

namespace TestEstimatorMergeProbs
{
	struct testPoissonMergeRejections;
}

namespace Estimation
{
namespace Merge
{
	namespace BarcodesParsing
	{
		class BarcodesParser;
	}

	class RealBarcodesMergeStrategy : public MergeStrategyBase
	{
		friend struct TestEstimator::testBarcodesFile;
		friend struct TestEstimator::testUmigsIntersection;
		friend struct TestEstimator::testFillDistances;
		friend struct TestEstimator::testRealNeighboursCbs;
		friend struct TestEstimator::testRealNeighbours;
		friend struct TestEstimatorMergeProbs::testPoissonMergeRejections;

	private:
		std::shared_ptr<BarcodesParsing::BarcodesParser> _barcodes_parser;

	private:
		ul_list_t get_real_neighbour_cbs(const CellsDataContainer &container, size_t base_cell_ind) const;

	protected:
		long get_merge_target(const CellsDataContainer &container, size_t base_cell_ind) const;
		virtual long get_best_merge_target(const CellsDataContainer &container, size_t base_cell_ind, const ul_list_t &neighbour_cells) const;
		virtual unsigned get_max_merge_dist(unsigned min_real_cb_dist) const;

	public:
		RealBarcodesMergeStrategy(const std::string &barcodes_filename, const boost::property_tree::ptree &config);

		virtual std::string merge_type() const override;
	};
}


}
