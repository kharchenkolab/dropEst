#pragma once

#include "SpacerFinder.h"
#include "TagsFinderBase.h"

namespace TestTagsSearch
{
	struct test1;
	struct test2;
	struct test3;
}

namespace TagsSearch
{
	class IndropV1TagsFinder : public TagsFinderBase
	{
		friend struct TestTagsSearch::test1;
		friend struct TestTagsSearch::test2;
		friend struct TestTagsSearch::test3;
	private:
		SpacerFinder _spacer_finder;

	protected:
		virtual Tools::ReadParameters parse(const std::string &r1_seq, const std::string &r1_quality,
		                                    const SpacerFinder::spacer_pos_t &spacer_pos);

		bool parse_fastq_record(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params) override;
		std::string get_additional_stat(long total_reads_read) const override;

	public:
		IndropV1TagsFinder(const std::vector<std::string> &fastq_filenames,
		                   const boost::property_tree::ptree &spacer_config, const boost::property_tree::ptree &config,
		                   const std::shared_ptr<ConcurrentGzWriter> &writer, bool save_stats, bool save_read_params);
	};
}

