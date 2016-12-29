#include "MergeStrategyFactory.h"

#include "BrokenRealBarcodesMergeStrategy.h"
#include "DummyMergeStrategy.h"
#include "PoissonRealBarcodesMergeStrategy.h"
#include "RealBarcodesMergeStrategy.h"
#include "SimpleMergeStrategy.h"
#include "PoissonSimpleMergeStrategy.h"

namespace Estimation
{
namespace Merge
{
	std::shared_ptr<MergeStrategyAbstract>
	MergeStrategyFactory::get(const boost::property_tree::ptree &config, bool merge_tags, const std::string &barcodes_filename)
	{
		if (!merge_tags)
			return std::shared_ptr<MergeStrategyAbstract>(new DummyMergeStrategy(config));

		std::string merge_type(config.get<std::string>("merge_type", ""));
		if (barcodes_filename == "")
		{
			if (merge_type == "poisson")
				return std::shared_ptr<MergeStrategyAbstract>(new PoissonSimpleMergeStrategy(config));

			return std::shared_ptr<MergeStrategyAbstract>(new SimpleMergeStrategy(config));
		}

		if (merge_type == "")
			return std::shared_ptr<MergeStrategyAbstract>(new RealBarcodesMergeStrategy(barcodes_filename, config));

		if (merge_type == "poisson")
			return std::shared_ptr<MergeStrategyAbstract>(new PoissonRealBarcodesMergeStrategy(barcodes_filename, config));

		if (merge_type == "broken")
			return std::shared_ptr<MergeStrategyAbstract>(new BrokenRealBarcodesMergeStrategy(barcodes_filename, config));

		throw std::runtime_error("Unknown merge type: " + merge_type);
	}
}
}