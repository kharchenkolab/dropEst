#pragma once

#include "TagsFinderBase.h"

namespace TagsSearch
{
	class IClipTagsFinder : public TagsFinderBase
	{
	private:
		const size_t _barcode_length;
		const size_t _umi_length;
		size_t _cant_parse_num;

	protected:
		bool parse_fastq_records(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params) override;

		std::string get_additional_stat(long total_reads_read) const override;

		std::string parse_cb(const std::string &sequence) const;
		std::string parse_umi(const std::string &sequence) const;
		std::string trim_barcodes(const std::string &sequence) const;

	public:
		IClipTagsFinder(const std::vector<std::string> &fastq_filenames,
		                const boost::property_tree::ptree &barcodes_config,
		                const boost::property_tree::ptree &config,
		                const std::shared_ptr<ConcurrentGzWriter> &writer, bool save_stats, bool save_read_params);
	};
}

