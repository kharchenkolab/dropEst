#pragma once

#include "SpacerFinder.h"
#include "TagsFinderBase.h"


namespace TagsSearch
{
	class SplitSeqTagsFinder : public TagsFinderBase
	{
	private:
		const std::vector<size_t> _barcode_starts;
		const std::vector<size_t> _barcode_lengths;
		const size_t _umi_start;
		const size_t _umi_length;
		size_t _min_barcode_read_length;
		size_t _short_seq_read_num;

	protected:
		bool parse_fastq_records(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params) override;
		std::string get_additional_stat(long total_reads_read) const override;

		std::string parse_cell_barcode(const std::string &sequence);
		std::string parse_umi_barcode(const std::string &sequence);

	public:
		SplitSeqTagsFinder(const std::vector<std::string> &fastq_filenames,
		                   const boost::property_tree::ptree &barcode_config,
		                   const boost::property_tree::ptree &processing_config,
		                   const std::shared_ptr<ConcurrentGzWriter> &writer,
		                   bool save_stats, bool save_read_params);
	};
}

