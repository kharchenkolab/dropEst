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
			const std::string cb;
			const std::string umi;
			const std::string gene;

			BamTags();
			BamTags(const boost::property_tree::ptree &config);
		};
	}
}