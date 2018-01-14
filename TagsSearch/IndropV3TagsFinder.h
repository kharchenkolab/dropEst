#pragma once

#include <TagsSearch/Counters/TwoBarcodesCounter.h>
#include "TagsFinderBase.h"

namespace TagsSearch
{
	class IndropV3TagsFinder : public TagsFinderBase
	{
	private:
		using len_t = std::string::size_type;

		const len_t barcode1_length;
		const len_t barcode2_length;
		const len_t umi_length;
		const len_t trim_tail_length;

		TwoBarcodesCounter _counter;


	private:
		std::string parse_cb(const std::string &cb1_seq, const std::string &cb2_seq) const;

	protected:
		bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) override;

		std::string get_additional_stat(long total_reads_read) const override;

	public:
		IndropV3TagsFinder(const std::vector<std::string> &fastq_filenames,
		                   const boost::property_tree::ptree &barcodes_config,
		                   const boost::property_tree::ptree &processing_config,
		                   const std::shared_ptr<ConcurrentGzWriter> &writer,
		                   bool save_stats, bool save_read_params);

		std::string parse_umi(const std::string &cb2_seq) const;
	};
}