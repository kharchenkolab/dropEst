#include "FastQReader.h"

#include <sstream>
#include <stdexcept>

#include <Tools/Logs.h>

namespace TagsSearch
{
	FastQReader::FastQReader(const std::string &filename)
		: _filename(filename)
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

	FastQReader::FastQRecord FastQReader::get_next_record()
	{
		FastQRecord record;
		if (!std::getline(this->_in_fstream, record.id).good())
			return record;

		if (record.id.at(0) != '@')
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq malformed!");

		if (!std::getline(this->_in_fstream, record.sequence).good())
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq ended prematurely!");

		if (!std::getline(this->_in_fstream, record.description).good())
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq ended prematurely!");

		if (!std::getline(this->_in_fstream, record.quality).good())
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': fastq ended prematurely!");

		if (record.sequence.length() != record.quality.length())
			throw std::runtime_error("File '" + this->_filename + "', read '" + record.id + "': different lengths of the sequence and the quality string!");

		return record;
	}

	const std::string &FastQReader::filename() const
	{
		return this->_filename;
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
