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
			, read_type(config.get<std::string>("BamTags.Type.tag", ""))
			, intronic_read_value(config.get<std::string>("BamTags.Type.intronic", ""))
			, intergenic_read_value(config.get<std::string>("BamTags.Type.intergenic", ""))
		{
			if (!this->read_type.empty() && this->intronic_read_value.empty())
				throw std::runtime_error("You have to specify tag values to be able to parse info about read types (see conf_desc.xml \"Estimation/BamTags/Type/\")");
		}

		BamTags::BamTags()
			: cell_barcode("")
			, umi("")
			, gene("")
			, cb_quality("")
			, umi_quality("")
		{}
	}
}