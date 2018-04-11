#pragma once

#include <string>
#include <boost/property_tree/ptree.hpp>

namespace Estimation
{
	namespace BamProcessing
	{
		class BamTags
		{
		public:
			const std::string cell_barcode;
			const std::string cell_barcode_raw;
			const std::string umi;
			const std::string umi_raw;
			const std::string gene;
			const std::string cell_barcode_quality;
			const std::string umi_quality;

			const std::string read_type;
			const std::string intronic_read_value;
			const std::string intergenic_read_value;
			const std::string exonic_read_value;

			const std::string intronic_read_value_out;
			const std::string intergenic_read_value_out;
			const std::string exonic_read_value_out;

			BamTags() = default;

			explicit BamTags(const boost::property_tree::ptree &config);
		};
	}
}