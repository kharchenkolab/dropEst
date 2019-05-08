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
	struct testTrimming;
	struct testValidation;
}

namespace Tools
{
	class ReadParameters;
}

namespace TagsSearch
{
	/// Main class, which should be extended for each new protocol
	/// Performs parallel reading of set of FastQ files to extract information about Cell Barcode, UMI and gene sequence for each read
	class TagsFinderBase
	{
		friend struct TestTagsSearch::test1;
		friend struct TestTagsSearch::testTrimming;
		friend struct TestTagsSearch::testValidation;

	public:
		using s_counter_t = std::unordered_map<std::string, int>;

	protected:
		using len_t = std::string::size_type;

	private:
		const bool _save_stats;
		const bool _save_read_params;
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
		const char _barcode_phred_threshold;
		const char _trim_phred_threshold;
		const char _gene_phred_threshold;
		const int _leading_trim_length;
		const int _trailing_trim_length;
		const double _max_g_fraction;
		const std::string _poly_a;
		const Tools::ReverseComplement _rc;

		TrimsCounter _trims_counter;

	private:
		static std::string get_file_uid(long random_seed = -1);

		/// Parse next block of fastq records from the readers and compile them to a single output fastq records with its parameters.
		///
		/// \param gene_record (output) compiled record with gene read
		/// \param params (output) read parameters for the compiled record
		/// \return false in the next cases: (1) fastq files are ended, (2) record doesn't pass quality threshold,
		/// (3) can't parse barcodes from the record; true otherwise
		bool get_next_output_record(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &params);

		/// Parse next bunch of fastq records and its read parameters, compile them to string output and
		/// place it to the corresponding queues
		///
		/// \param number_of_iterations number of reading iterations, each of them read records_bunch_size records,
		/// compile it to a single string and put as a single record to queues
		/// \param records_bunch_size number of records per iteration
		void parse_next_bunch(size_t number_of_iterations = 10000, size_t records_bunch_size = 5000);

		/// Run processing of fastq file for one thread
		void run_thread();

		bool validate(const FastQReader::FastQRecord& record) const;
		bool trim(FastQReader::FastQRecord& record) const;

	protected:

		/// Get next block from each of fastq files, parse read parameters from them and generate single new fastq
		/// record which will be saved for further transcriptome alignment.
		/// Must be overwritten for each new protocol.
		///
		/// \param record (output) fastq record with gene read
		/// \param read_params (output) parsed read parameters
		///
		/// \return false if end of fastq files is reached, true otherwise
		virtual bool parse_fastq_records(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params) = 0;

		/// Trim poly-A tail from the gene read
		///
		/// \param sequence gene read sequence, which should be trimmed (inplace)
		/// \param quality gene read quality sequence, which should be trimmed (inplace)
		/// \param barcode_tail tail of cell barcode. If parameter is provided, pipeline tries to find reverse-complete
		/// tail sequence for trimming
		void trim_poly_a(std::string &sequence, std::string &quality, const std::string &barcode_tail="");

		/// Convert some additional statistics, collected during parsing, to print them in logfile.
		/// Can be safely ignored during implementation of new protocol parsers.
		///
		/// \param total_reads_read total number of reads in the file
		virtual std::string get_additional_stat(long total_reads_read) const;

		/// Get reader for the specific fastq file
		///
		/// \param index fastq file index
		/// \return corresponding fastq reader
		FastQReader& fastq_reader(size_t index);

	public:
		TagsFinderBase(const std::vector<std::string> &fastq_filenames,
		               const boost::property_tree::ptree &processing_config, const std::shared_ptr<ConcurrentGzWriter> &writer,
		               bool save_stats, bool save_read_params);

		/// Create workers, which perform parallel processing of the fastq files
		/// \param number_of_threads number of workers
		void run(int number_of_threads);
		const s_counter_t& num_reads_per_cb() const;
		std::string results_to_string() const;
	};
}