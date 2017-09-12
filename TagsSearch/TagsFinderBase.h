#pragma once

#include "Counters/TrimsCounter.h"
#include "FastQReader.h"
#include "SpacerFinder.h"
#include "Tools/UtilFunctions.h"

#include <string>

#include <boost/property_tree/ptree.hpp>
#include <unordered_map>

namespace TestTagsSearch
{
	struct test1;
}

namespace Tools
{
	class ReadParameters;
}

namespace TagsSearch
{
	class TagsFinderBase
	{
		friend struct TestTagsSearch::test1;
	public:
		typedef std::unordered_map<std::string, int> s_counter_t;

	protected:
		typedef std::string::size_type len_t;

	private:
		const bool _save_stats;
		const std::string _file_uid;

		s_counter_t _num_reads_per_cb;
		long _total_reads_read;
		long _parsed_reads;
		bool _file_ended;

	protected:
		const size_t _max_reads;
		const unsigned _min_read_len;
		const std::string poly_a;
		const Tools::ReverseComplement rc;

		TrimsCounter _trims_counter;

	public:
		long long read_time, parse_time, construct_time;

	private:
		static std::string get_file_uid(long random_seed);

	protected:
		void trim(const std::string &barcodes_tail, std::string &sequence, std::string &quality);
		virtual bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) = 0;
		virtual std::string get_additional_stat(long total_reads_read) const = 0;

	public:
		TagsFinderBase(const boost::property_tree::ptree &processing_config, bool save_stats);

		virtual bool get_next_record(FastQReader::FastQRecord& record);
		const s_counter_t& num_reads_per_cb() const;
		bool file_ended() const;
		std::string results_to_string() const;
	};
}