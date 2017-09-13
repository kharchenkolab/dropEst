#pragma once

#include "IndropV3TagsFinder.h"

namespace TagsSearch
{
	class IndropV3LibsTagsFinder : public IndropV3TagsFinder
	{
	private:
		const std::string library_tag;
		const unsigned max_lib_tag_ed;

		FastQReader _lib_tag_reader;

	protected:
		virtual bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) override;

	public:
		IndropV3LibsTagsFinder(const std::string &barcode1_fastq_name, const std::string &barcode2_fastq_name,
		                       const std::string &gene_fastq_name, const std::string &library_fastq_name,
		                       const std::string &library_tag, const boost::property_tree::ptree &barcodes_config,
		                       const boost::property_tree::ptree &config, TextWriter &&writer, bool save_stats);
	};
}
