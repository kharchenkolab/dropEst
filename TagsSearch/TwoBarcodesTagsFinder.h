#pragma once

#include <TagsSearch/Counters/TwoBarcodesCounter.h>
#include "TagsFinderBase.h"

namespace TagsSearch
{
	class TwoBarcodesTagsFinder : public TagsFinderBase
	{
	private:
		typedef std::string::size_type len_t;

		const len_t barcode1_length;
		const len_t barcode2_length;
		const len_t umi_length;
		const len_t trim_tail_length;

		TwoBarcodesCounter counter;

	public:
		TwoBarcodesTagsFinder(const std::shared_ptr<FilesProcessor> &files_processor,
							  const boost::property_tree::ptree &barcodes_config,
							  const boost::property_tree::ptree &config);

	protected:
		virtual bool parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params) override;
		virtual std::string get_additional_stat(long total_reads_read) const override;
	};
}