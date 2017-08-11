#pragma once

#include "Counters/TrimsCounter.h"
#include "FilesProcessor.h"
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
	class FilesProcessor;

	class TagsFinderBase
	{
		friend struct TestTagsSearch::test1;
	public:
		typedef std::unordered_map<std::string, int> s_counter_t;

	protected:
		typedef std::string::size_type len_t;

	private:
		s_counter_t _num_reads_per_cb;

	protected:
		const size_t max_reads;
		const unsigned min_read_len;
		const std::string poly_a;
		const Tools::ReverseComplement rc;

		const std::shared_ptr<FilesProcessor> _files_processor;
		TrimsCounter _trims_counter;

	private:
		std::string results_to_string(long total_reads_read) const;

	protected:
		void trim(const std::string &barcodes_tail, std::string &sequence, std::string &quality);
		virtual bool parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params) = 0;
		virtual std::string get_additional_stat(long total_reads_read) const = 0;

	public:
		TagsFinderBase(std::shared_ptr<FilesProcessor> files_processor, const boost::property_tree::ptree &processing_config);

		virtual void run(bool save_reads_names, bool save_stats);
		const s_counter_t& num_reads_per_cb() const;
	};
}