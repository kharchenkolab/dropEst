#pragma once

#include "MergeStrategyAbstract.h"

namespace Estimation
{
namespace Merge
{
	class MergeStrategyFactory
	{
	public:
		static std::shared_ptr<MergeStrategyAbstract> get(const boost::property_tree::ptree &config,
														  bool merge_tags, const std::string &barcodes_filename);
	};
}
}