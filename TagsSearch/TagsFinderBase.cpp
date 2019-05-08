#include "TagsFinderBase.h"

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"

#include <thread>

using Tools::ReadParameters;

namespace TagsSearch
{
	TagsFinderBase::TagsFinderBase(const std::vector<std::string> &fastq_filenames,
	                               const boost::property_tree::ptree &processing_config,
	                               const std::shared_ptr<ConcurrentGzWriter> &writer,
	                               bool save_stats, bool save_read_params)
		: _save_stats(save_stats)
		, _save_read_params(save_read_params)
		, _barcode_phred_threshold(ReadParameters::quality_to_phred(processing_config.get<int>("min_barcode_quality", 0)))
		, _trim_phred_threshold(ReadParameters::quality_to_phred(processing_config.get<int>("trim_quality", 0)))
		, _gene_phred_threshold(ReadParameters::quality_to_phred(processing_config.get<int>("min_median_quality", 0)))
		, _leading_trim_length(processing_config.get<int>("leading_trim", 0))
		, _trailing_trim_length(processing_config.get<int>("trailing_trim", 0))
		, _max_g_fraction(processing_config.get<double>("max_g_fraction", 1.0))
		, _file_uid(TagsFinderBase::get_file_uid())
		, _total_reads_read(0)
		, _low_quality_reads(0)
		, _parsed_reads(0)
		, _file_ended(false)
		, _reading_in_progress(false)
		, _min_read_len(processing_config.get<unsigned>("min_align_length", 10))
		, _poly_a(processing_config.get<std::string>("poly_a_tail", "AAAAAAAA"))
		, _trims_counter()
		, _fastq_writer(writer)
	{
		for (auto &&filename : fastq_filenames)
		{
			this->_fastq_readers.emplace_back(std::make_shared<FastQReader>(filename));
		}

		if (save_read_params)
		{
			this->_params_writer = std::make_shared<ConcurrentGzWriter>(writer->base_filename(), "params.gz", 0);
		}
	}

	bool TagsFinderBase::get_next_output_record(FastQReader::FastQRecord &gene_record, ReadParameters &params)
	{
		if (this->_file_ended)
			return false;

		if (!this->parse_fastq_records(gene_record, params))
		{
			this->_file_ended = true;
			return false;
		}

		if (this->_total_reads_read % 5000000 == 0 && this->_total_reads_read > 0)
		{
			L_TRACE << "Total " << this->_total_reads_read << " read (" << this->_parsed_reads << " parsed, "
			        << (this->_parsed_reads - this->_low_quality_reads) << " passed quality threshold)";
		}
		this->_total_reads_read++;

		if (params.is_empty() || gene_record.sequence.length() < this->_min_read_len)
			return false;

		++this->_parsed_reads;

		if (!params.pass_quality_threshold() || !this->validate(gene_record) || !this->trim(gene_record))
		{
			this->_low_quality_reads++;
			return false;
		}

		std::string read_prefix = "@" + this->_file_uid + std::to_string(this->_total_reads_read);

		gene_record.id = this->_save_read_params ? read_prefix : params.encoded_id(read_prefix);

		if (this->_save_stats)
		{
			this->_num_reads_per_cb[params.cell_barcode()]++;
		}

		return true;
	}

	std::string TagsFinderBase::results_to_string() const
	{
		std::stringstream ss;
		ss << " (" << this->_total_reads_read << " reads)\n"
		   << this->get_additional_stat(this->_total_reads_read) << "\n"
		   << this->_trims_counter.print();

		return ss.str();
	}

	void TagsFinderBase::trim_poly_a(std::string &sequence, std::string &quality, const std::string &barcode_tail)
	{
		if (sequence.length() != quality.length())
			throw std::runtime_error("Read has different lengths of sequence and quality string: '" +
											 sequence + "', '" + quality + "'");

		len_t trim_pos = sequence.length();
		if (barcode_tail.length() > 0)
		{
			// attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
			// RC of UMI+second barcode (up to a length r1_rc_length - spacer_finder parameter)
			std::string rcb = this->_rc.rc(barcode_tail);

			len_t rc_pos = sequence.find(rcb);
			if (rc_pos != std::string::npos)
			{
				trim_pos = rc_pos;
				this->_trims_counter.inc(TrimsCounter::RC);
			}
		}

		if (trim_pos == sequence.length())
		{
			// attempt 2: find polyA block
			len_t rc_pos = sequence.find(this->_poly_a);
			if (rc_pos != std::string::npos)
			{
				trim_pos = rc_pos;
				this->_trims_counter.inc(TrimsCounter::POLY_A);
			}
		}

		// attempt 3: trim trailing As
		bool a_trim = false;
		len_t skip_count = 0;
		while (trim_pos > 0 && (sequence.at(trim_pos - 1) == 'A' || sequence.at(trim_pos - 1) == 'N'))
		{
			trim_pos--;
			skip_count++;
			a_trim = true;
		}
		if (a_trim)
		{
			this->_trims_counter.inc(TrimsCounter::A_TRIM);
		}

		//attempt 4: apply
		if (sequence.length() != trim_pos)
		{
			sequence = sequence.substr(0, trim_pos);
			quality = quality.substr(0, trim_pos);
		}
		else
		{
			this->_trims_counter.inc(TrimsCounter::NO_TRIM);
		}
	}

	std::string TagsFinderBase::get_additional_stat(long total_reads_read) const
	{
		return "";
	}

	const TagsFinderBase::s_counter_t& TagsFinderBase::num_reads_per_cb() const
	{
		return this->_num_reads_per_cb;
	}

	std::string TagsFinderBase::get_file_uid(long random_seed)
	{
		if (random_seed < 0)
		{
			random_seed = time(nullptr);
		}

		srand(unsigned(random_seed));

		const int range = 'Z' - 'A';
		std::string res;
		for (int i = 0; i < 4; ++i)
		{
			res.append(1, char(rand() % range + 'A'));
		}
		return res;
	}

	void TagsFinderBase::parse_next_bunch(size_t number_of_iterations, size_t records_bunch_size)
	{
		for (size_t i = 0; i < number_of_iterations; ++i)
		{
			if (this->_file_ended)
				break;

			std::string records_bunch, params_bunch;
			unsigned record_id;
			for (record_id = 0; record_id < records_bunch_size; ++record_id)
			{
				if (this->_file_ended)
					break;

				FastQReader::FastQRecord record;
				ReadParameters params;
				if (!this->get_next_output_record(record, params))
				{
					record_id--;
					continue;
				}

				records_bunch += record.to_string();
				params_bunch += params.to_string(record.id) + "\n";
			}

			if (!records_bunch.empty())
			{
				this->_fastq_writer->enqueue_lines(records_bunch, record_id);

				if (this->_save_read_params)
				{
					this->_params_writer->enqueue_lines(params_bunch, record_id);
				}
			}
		}
	}

	void TagsFinderBase::run_thread()
	{
		while (true)
		{
			if (this->_file_ended && this->_fastq_writer->empty())
				break;

			// Part 1. Multithreaded over files. Load some content from each of fastq files to memory.
			for (auto &reader : this->_fastq_readers)
			{
				reader->try_read_records_to_cash();
			}

			// Part 2. Single thread. Get cached records from fasq readers and parse them to gene records and
			// (optionally) read parameters, caching results to file writers
			if (!this->_file_ended)
			{
				if (!this->_reading_in_progress.exchange(true))
				{
					this->parse_next_bunch();
					if (this->_file_ended)
					{
						L_TRACE << this->results_to_string();
						Tools::trace_time("Reading completed");
						L_TRACE << "Writing the rest lines";
					}

					this->_reading_in_progress = false;
				}
			}

			// Part 3.1. Multithreaded. Gzip cached lines with gene records
			this->_fastq_writer->gzip_cached_lines(this->_file_ended);

			// Part 3.2. Multithreaded. Gzip cached lines with read parameters
			if (this->_save_read_params)
			{
				this->_params_writer->gzip_cached_lines(_file_ended);
			}

			// Part 4.1. Single thread. Write all cached gzipped gene records.
			this->_fastq_writer->flush_writing_queue();

			// Part 4.2. Single thread. Write all cached gzipped read parameters.
			if (this->_save_read_params)
			{
				this->_params_writer->flush_writing_queue();
			}
		}
	}

	void TagsFinderBase::run(int number_of_threads)
	{
		std::vector<std::thread> tasks;
		for (int thread_num = 0; thread_num < number_of_threads; ++thread_num)
		{
			tasks.emplace_back([this]{this->run_thread();});
		}

		for (auto &task : tasks)
		{
			task.join();
		}
	}

	FastQReader &TagsFinderBase::fastq_reader(size_t index)
	{
		return *this->_fastq_readers.at(index);
	}

	bool TagsFinderBase::validate(const FastQReader::FastQRecord &record) const
	{
		if (this->_gene_phred_threshold <= ReadParameters::quality_offset)
			return true;

		double n_low = 0;
		for (char qual : record.quality)
		{
			n_low += (qual < this->_gene_phred_threshold);
		}

		if (n_low / record.quality.size() > 0.5)
			return false;

		double n_g = 0;
		for (char nuc : record.sequence)
		{
			n_g += (nuc == 'G' || nuc == 'N');
		}

		return (n_g / record.quality.size() < this->_max_g_fraction);
	}

	bool TagsFinderBase::trim(FastQReader::FastQRecord &record) const
	{
		if (this->_trim_phred_threshold <= ReadParameters::quality_offset)
			return true;

		int trim_start = -1;
		for (int i = 0; i < std::min(int(record.sequence.size()), this->_leading_trim_length); ++i)
		{
			if (record.quality.at(static_cast<unsigned>(i)) < this->_trim_phred_threshold)
			{
				trim_start = i;
			}
		}

		trim_start++;
		long trim_end = static_cast<int>(record.sequence.size());
		for (long i = record.sequence.size() - 1; i >= std::max(long(record.sequence.size() - this->_trailing_trim_length), 0L); --i)
		{
			if (record.quality.at(static_cast<unsigned long>(i)) < this->_trim_phred_threshold)
			{
				trim_end = i;
			}
		}

		auto trimmed_length = static_cast<size_t>(trim_end - trim_start);
		if (trimmed_length < this->_min_read_len)
			return false;

		record.quality = record.quality.substr(static_cast<size_t>(trim_start), trimmed_length);
		record.sequence = record.sequence.substr(static_cast<size_t>(trim_start), trimmed_length);
		return true;
	}
}
