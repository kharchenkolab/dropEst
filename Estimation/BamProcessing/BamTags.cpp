#include "BamTags.h"

namespace Estimation
{
	namespace BamProcessing
	{
		BamTags::BamTags(const boost::property_tree::ptree &config)
			: cell_barcode(config.get<std::string>("BamTags.cb", "CB"))
			, cell_barcode_raw(config.get<std::string>("BamTags.cb_raw", "CR"))
			, umi(config.get<std::string>("BamTags.umi", "UB"))
			, umi_raw(config.get<std::string>("BamTags.umi_raw", "UR"))
			, gene(config.get<std::string>("BamTags.gene", "GX"))
			, cell_barcode_quality(config.get<std::string>("BamTags.cb_quality", "CQ"))
			, umi_quality(config.get<std::string>("BamTags.umi_quality", "UQ"))
			, read_type(config.get<std::string>("BamTags.Type.tag", ""))
			, intronic_read_value(config.get<std::string>("BamTags.Type.intronic", ""))
			, intergenic_read_value(config.get<std::string>("BamTags.Type.intergenic", ""))
			, exonic_read_value(config.get<std::string>("BamTags.Type.exonic", ""))
			, intronic_read_value_out(config.get<std::string>("BamTags.Type.intronic", "INTRONIC"))
			, intergenic_read_value_out(config.get<std::string>("BamTags.Type.intergenic", "INTERGENIC"))
			, exonic_read_value_out(config.get<std::string>("BamTags.Type.exonic", "EXONIC"))
		{
			if (!this->read_type.empty() && this->intronic_read_value.empty())
				throw std::runtime_error("You have to specify tag values to be able to parse info about read types (see conf_desc.xml \"Estimation/BamTags/Type/\")");
		}
	}
}