#include "BamTags.h"

namespace Estimation
{
	namespace BamProcessing
	{
		BamTags::BamTags(const boost::property_tree::ptree &config)
			: cell_barcode(config.get<std::string>("BamTags.cb", "CB"))
			, umi(config.get<std::string>("BamTags.umi", "UB"))
			, gene(config.get<std::string>("BamTags.gene", "GX"))
			, cb_quality(config.get<std::string>("BamTags.cb_quality", "CQ"))
			, umi_quality(config.get<std::string>("BamTags.umi_quality", "UQ"))
		{}

		BamTags::BamTags()
			: cell_barcode("")
			, umi("")
			, gene("")
			, cb_quality("")
			, umi_quality("")
		{}
	}
}