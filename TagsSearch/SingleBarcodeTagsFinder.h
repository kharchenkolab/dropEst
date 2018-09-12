#pragma once

#include "TagsFinderBase.h"

namespace TagsSearch
{
	class SingleBarcodeTagsFinder : public TagsFinderBase
	{
	private:
		using len_t = std::string::size_type;

		const len_t barcode_length;
		const len_t umi_length;
		const len_t trim_tail_length;

	private:
		std::string parse_cb(const std::string &cb_seq) const;

	protected:
		bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) override;

	public:
		SingleBarcodeTagsFinder(const std::vector<std::string> &fastq_filenames,
		                        const boost::property_tree::ptree &barcodes_config,
		                        const boost::property_tree::ptree &processing_config,
		                        const std::shared_ptr<ConcurrentGzWriter> &writer,
		                        bool save_stats, bool save_read_params);

		std::string parse_umi(const std::string &cb_seq) const;
	};
}