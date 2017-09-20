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
			const std::string umi;
			const std::string gene;
			const std::string cb_quality;
			const std::string umi_quality;

			BamTags();
			BamTags(const boost::property_tree::ptree &config);
		};
	}
}