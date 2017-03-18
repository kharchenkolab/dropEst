#pragma once

#include "IndropV3TagsFinder.h"

namespace TagsSearch
{
	class IndropV3LibsTagsFinder : public IndropV3TagsFinder
	{
	private:
		std::string library_tag;
		const unsigned max_lib_tag_ed;

	protected:
		virtual bool parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params) override;

	public:
		IndropV3LibsTagsFinder(const std::shared_ptr<FilesProcessor> &files_processor, const std::string &library_tag,
		                       unsigned max_lib_tag_ed, const boost::property_tree::ptree &barcodes_config,
		                       const boost::property_tree::ptree &config);
	};
}
