#pragma once

#include "SpacerFinder.h"
#include "Counters/TrimsCounter.h"

#include <string>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/property_tree/ptree.hpp>

namespace TestTagsSearch
{
	struct test1;
	struct test2;
	struct test3;
}

namespace Tools
{
	class ReadParameters;
}

namespace TagsSearch
{
	class FilesProcessor;

	class TagsFinder
	{
		friend struct TestTagsSearch::test1;
		friend struct TestTagsSearch::test2;
		friend struct TestTagsSearch::test3;

	private:
		typedef std::string::size_type len_t;

		long max_reads;
		unsigned min_align_len;
		std::string poly_a;

		SpacerFinder spacer_finder;
		TrimsCounter trims_counter;

	private:
		static bool read_blocks(FilesProcessor &files_processor, long total_reads_read,
								std::string &r1_seq, std::string &r2_id, std::string &r2_seq,
								std::string &r2_description, std::string &r2_quality);

		std::string results_to_string(long total_reads_read) const;
		len_t get_trim_position(len_t spacer_end, const std::string &r1_seq, const std::string &r2_seq);
		Tools::ReadParameters fill_parameters(const std::string &r1_seq, const std::string &r2_id,
											  std::string &r2_seq, const std::string &r2_description,
											  std::string &r2_quality_str);

	public:
		TagsFinder()
		{}

		TagsFinder(const SpacerFinder &spacer_finder, const boost::property_tree::ptree &config);

		void run(const std::string &r1_filename, const std::string &r2_filename, const std::string &base_name,
				 bool save_reads_names);
	};
}