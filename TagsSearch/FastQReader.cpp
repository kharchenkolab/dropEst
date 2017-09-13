#include "FastQReader.h"

#include <sstream>
#include <stdexcept>

#include <Tools/Logs.h>

namespace TagsSearch
{
	FastQReader::FastQReader(const std::string &filename)
		: _filename(filename)
		, _file_ended(false)
		, _records(10000)
	{
		if (filename.empty())
			return;

		this->_in_file.open(filename, std::ios_base::in | std::ios_base::binary);
		if (!this->_in_file)
			throw std::runtime_error("Can't open fastq file '" + filename + "'");

		if (boost::ends_with(filename, ".gz") || boost::ends_with(filename, ".gzip"))
		{
			this->_in_fstream.push(boost::iostreams::gzip_decompressor());
		}
		this->_in_fstream.push(this->_in_file);
	}

//	void FastQReader::write_read_params(const std::string &id_prefix, const Tools::ReadParameters &read_params)
//	{
//		this->out_reads_zip << read_params.encoded_params(id_prefix) << "\n";
//	}

	bool FastQReader::get_next_record_unsafe(FastQRecord &record)
	{
		if (!this->get_next_line_unsafe(record.id))
			return false;

		if (record.id.at(0) != '@')
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq malformed!");

		if (!this->get_next_line_unsafe(record.sequence))
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq ended prematurely!");

		if (!this->get_next_line_unsafe(record.description))
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq ended prematurely!");

		if (!this->get_next_line_unsafe(record.quality))
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq ended prematurely!");

		if (record.sequence.length() != record.quality.length())
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': different lengths of the sequence and the quality string!");

		return true;
	}

	const std::string &FastQReader::filename() const
	{
		return this->_filename;
	}

	bool FastQReader::get_next_line_unsafe(std::string &line)
	{
		bool success = std::getline(this->_in_fstream, line).good();
		this->_file_ended = !success;
		return success;
	}

	void FastQReader::try_read_records_to_cash()
	{
		if (this->_file_ended)
			return;

		std::unique_lock<mutex_t> lock(this->_read_mutex, std::try_to_lock);
		if (!lock)
			return;

		while (!this->_records.full())
		{
			FastQRecord record;
			if (!this->get_next_record_unsafe(record))
				break;

			this->_records.push(record);
		}
	}

	bool FastQReader::get_next_record(FastQReader::FastQRecord &record)
	{
		while (true)
		{
			if (this->_records.pop(record))
				break;

			if (this->_file_ended && this->_records.empty())
				return false;

			this->try_read_records_to_cash();
		}

		return true;
	}

	FastQReader::FastQRecord::FastQRecord(const std::string &id, const std::string &sequence,
	                                      const std::string &description, const std::string &quality)
		: id(id)
		, sequence(sequence)
		, description(description)
		, quality(quality)
	{}

	std::string FastQReader::FastQRecord::to_string() const
	{
		return this->id + "\n" + this->sequence + "\n" + this->description + "\n" + this->quality + "\n";
	}
}
