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
	class SpacerTagsFinder : public TagsFinderBase
	{
		friend struct TestTagsSearch::test1;
		friend struct TestTagsSearch::test2;
		friend struct TestTagsSearch::test3;
	private:
		SpacerFinder spacer_finder;

	protected:
		virtual Tools::ReadParameters parse(const std::string &r1_seq, const std::string &r1_quality,
				                                    const SpacerFinder::spacer_pos_t &spacer_pos);

		virtual bool parse_fastq_record(FilesProcessor::FastQRecord &gene_record, Tools::ReadParameters &read_params) override;
		virtual std::string get_additional_stat(long total_reads_read) const override;

	public:
		SpacerTagsFinder(std::shared_ptr<FilesProcessor> files_processor,
						 const boost::property_tree::ptree &spacer_config, const boost::property_tree::ptree &config);
	};
}

