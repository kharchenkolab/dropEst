#pragma once

#include <TagsSearch/Counters/TwoBarcodesCounter.h>
#include "TagsFinderBase.h"

namespace TagsSearch
{
	class IndropV3TagsFinder : public TagsFinderBase
	{
	private:
		typedef std::string::size_type len_t;

		const len_t barcode1_length;
		const len_t barcode2_length;
		const len_t umi_length;
		const len_t trim_tail_length;

		TwoBarcodesCounter _counter;

		FastQReader _barcode1_reader;
		FastQReader _barcode2_reader;
		FastQReader _gene_reader;


	private:
		std::string parse_cb(const std::string &cb1_seq, const std::string &cb2_seq) const;

	protected:
		virtual bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) override;
		virtual std::string get_additional_stat(long total_reads_read) const override;
		void skip_records_row();

	public:
		IndropV3TagsFinder(const std::string &barcode1_fastq_name, const std::string &barcode2_fastq_name,
		                   const std::string &gene_fastq_name, const boost::property_tree::ptree &barcodes_config,
		                   const boost::property_tree::ptree &processing_config, bool save_stats);

		std::string parse_umi(const std::string &cb2_seq) const;
	};
}