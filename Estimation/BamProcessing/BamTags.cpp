#include "BamTags.h"

namespace Estimation
{
	namespace BamProcessing
	{
		BamTags::BamTags(const boost::property_tree::ptree &config)
			: cb(config.get<std::string>("BamTags.cb", "CB"))
			, umi(config.get<std::string>("BamTags.umi", "UB"))
			, gene(config.get<std::string>("BamTags.gene", "GX"))
		{}

		BamTags::BamTags()
			: cb("")
			, umi("")
			, gene("")
		{}
	}
}