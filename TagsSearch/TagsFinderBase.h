#pragma once

#include "Counters/TrimsCounter.h"
#include "FastQReader.h"
#include "SpacerFinder.h"
#include "Tools/UtilFunctions.h"
#include <Tools/ScSpConcurrentQueue.h>
#include <Tools/BlockingConcurrentQueue.h>
#include "ConcurrentGzWriter.h"

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
		const bool _save_read_params;
		const int _quality_threshold;
		const std::string _file_uid;

		std::atomic<long> _total_reads_read;
		std::atomic<long> _low_quality_reads;
		std::atomic<long> _parsed_reads;
		std::atomic<bool> _file_ended;

		std::atomic<bool> _reading_in_progress;

		const std::shared_ptr<ConcurrentGzWriter> _fastq_writer;
		std::shared_ptr<ConcurrentGzWriter> _params_writer;
		std::vector<std::shared_ptr<FastQReader>> _fastq_readers;
		s_counter_t _num_reads_per_cb;

	protected:
		const unsigned _min_read_len;
		const std::string poly_a;
		const Tools::ReverseComplement rc;

		TrimsCounter _trims_counter;

	private:
		static std::string get_file_uid(long random_seed = -1);

		bool get_next_record(FastQReader::FastQRecord& record, Tools::ReadParameters &params);

		void read_bunch(size_t number_of_iterations = 10000, size_t records_bunch_size = 5000);
		void run_thread();

	protected:
		virtual bool parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) = 0;
		void trim(const std::string &barcodes_tail, std::string &sequence, std::string &quality);
		virtual std::string get_additional_stat(long total_reads_read) const = 0;
		FastQReader& fastq_reader(size_t index);

	public:
		TagsFinderBase(const std::vector<std::string> &fastq_filenames,
		               const boost::property_tree::ptree &processing_config, const std::shared_ptr<ConcurrentGzWriter> &writer,
		               bool save_stats, bool save_read_params);

		void run(int number_of_threads);
		const s_counter_t& num_reads_per_cb() const;
		std::string results_to_string() const;

		bool check_quality(const Tools::ReadParameters &parameters);
	};
}