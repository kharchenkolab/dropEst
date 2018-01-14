#pragma once

#include "IndropV3TagsFinder.h"

namespace TagsSearch
{
	class IndropV3LibsTagsFinder : public IndropV3TagsFinder
	{
	private:
		const std::string library_tag;
		const unsigned max_lib_tag_ed;

	protected:
		bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) override;

	public:
		IndropV3LibsTagsFinder(const std::vector<std::string> &fastq_filenames,
		                       const std::string &library_tag, const boost::property_tree::ptree &barcodes_config,
		                       const boost::property_tree::ptree &config, const std::shared_ptr<ConcurrentGzWriter> &writer,
		                       bool save_stats, bool save_read_params);
	};
}
