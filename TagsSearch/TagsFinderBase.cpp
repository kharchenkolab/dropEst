#include "TagsFinderBase.h"

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"

#include <thread>
#include <chrono>
#define chr std::chrono

namespace TagsSearch
{
	TagsFinderBase::TagsFinderBase(const std::vector<std::string> &fastq_filenames,
	                               const boost::property_tree::ptree &processing_config, TextWriter &&writer,
	                               bool save_stats)
		: _save_stats(save_stats)
		, _file_uid(TagsFinderBase::get_file_uid(42)) // TODO: return time(nullptr)
		, _total_reads_read(0)
		, _parsed_reads(0)
		, _file_ended(false)
		, _read_in_progress(false)
		, _write_in_progress(false)
		, _records(10000)
		, _gzipped(10000)
		, _min_read_len(processing_config.get<unsigned>("min_align_length", 10))
		, poly_a(processing_config.get<std::string>("poly_a_tail", "AAAAAAAA"))
		, _trims_counter()
		, _writer(std::move(writer))
	{
		for (auto &&filename : fastq_filenames)
		{
			this->_fastq_readers.emplace_back(std::make_shared<FastQReader>(filename));
		}
	}

	bool TagsFinderBase::get_next_record(FastQReader::FastQRecord& record)
	{
		if (this->_file_ended)
			return false;

		Tools::ReadParameters params;
		if (!this->parse_fastq_record(record, params))
		{
			this->_file_ended = true;
			return false;
		}

		if (++this->_total_reads_read % 1000000 == 0)
		{
			L_TRACE << "Total " << this->_total_reads_read << " read (" << this->_parsed_reads << " parsed)";
		}

		if (params.is_empty() || record.sequence.length() < this->_min_read_len)
			return false;

		++this->_parsed_reads;

		std::string read_prefix = "@" + this->_file_uid + std::to_string(this->_total_reads_read + 1); // TODO: remove "+1"
		record.id = params.encoded_id(read_prefix); // TODO: add save_read_names parameter // record.id = read_prefix;

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

	void TagsFinderBase::trim(const std::string &barcodes_tail, std::string &sequence, std::string &quality)
	{
		if (sequence.length() != quality.length())
			throw std::runtime_error("Read has different lengths of sequence and quality string: '" +
											 sequence + "', '" + quality + "'");

		len_t trim_pos = sequence.length();
		// attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
		// RC of UMI+second barcode (up to a length r1_rc_length - spacer_finder parameter)
		std::string rcb = this->rc.rc(barcodes_tail);

		len_t rc_pos = sequence.find(rcb);
		if (rc_pos != std::string::npos)
		{
			trim_pos = rc_pos;
			this->_trims_counter.inc(TrimsCounter::RC);
		}
		else
		{
			// attempt 2: find polyA block
			rc_pos = sequence.find(this->poly_a);
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
		return std::to_string(rand()) + char(rand() % 25 + 'A');
	}

	void TagsFinderBase::read_bunch(size_t number_of_iterations, size_t records_bunch_size)
	{
		for (size_t i = 0; i < number_of_iterations; ++i)
		{
			if (this->_file_ended)
				break;

			std::string records_bunch;
			for (size_t record_id = 0; record_id < records_bunch_size; ++record_id)
			{
				FastQReader::FastQRecord record;
				if (this->_file_ended)
					break;

				if (!this->get_next_record(record))
				{
					record_id--;
					continue;
				}

				records_bunch += record.to_string();
			}

			if (!records_bunch.empty())
			{
				this->_records.push(records_bunch);
			}
		}
	}

	void TagsFinderBase::gzip_all()
	{
		while (this->_file_ended || !this->_gzipped.full())
		{
			std::string records_bunch;

			if (!this->_records.pop(records_bunch))
				break;

			this->_gzipped.push(TextWriter::gzip(records_bunch));
		}
	}

	void TagsFinderBase::write_all()
	{
		if (this->_write_in_progress.exchange(true))
			return;

		std::string to_write, cur_text;
		while (this->_gzipped.pop(cur_text))
		{
			to_write += cur_text;
		}
		if (!to_write.empty())
		{
			this->_writer.write(to_write);
		}

		this->_write_in_progress = false;
	}

	void TagsFinderBase::run_thread()
	{
		while (true)
		{
			if (this->_file_ended && this->_records.empty() && this->_gzipped.empty())
				break;

			// Part 0. Multithreaded.
			for (auto &reader : this->_fastq_readers)
			{
				reader->try_read_records_to_cash();
			}

			// Part 1. Single thread.
			if (!this->_file_ended && !this->_read_in_progress.exchange(true))
			{
				this->read_bunch();
				if (this->_file_ended)
				{
					L_TRACE << this->results_to_string();
					Tools::trace_time("Reading completed");
					L_TRACE << "Writing the rest lines";
				}

				this->_read_in_progress = false;
			}

			// Part 2. Multithreaded.
			this->gzip_all();

			// Part 3. Single thread.
			this->write_all();
		}
	}

	void TagsFinderBase::run(int number_of_threads)
	{
		std::vector<std::thread> tasks;
		for (int thread_num = 0; thread_num < number_of_threads; ++thread_num)
		{
			tasks.push_back(std::thread([this]{this->run_thread();}));
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
}
