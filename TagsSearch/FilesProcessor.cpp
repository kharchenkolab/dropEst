#include "FilesProcessor.h"

#include <sstream>
#include <stdexcept>

#include <Tools/Logs.h>

namespace TagsSearch
{
	FilesProcessor::FilesProcessor(const std::vector<std::string> &filenames, const std::string &base_name,
								   bool save_reads_names)
			: filenames(filenames)
			, base_name(base_name)
			, reads_file_name(base_name + ".reads.gz")
			, current_file_reads_written(0)
			, out_file_index(1)
	{
		for (const auto &filename : filenames)
		{
			this->line_providers.push_back(std::make_shared<Tools::LineProvider>(filename, 8192));
		}

		std::string out_file_name = this->get_out_filename();
		this->out_file.open(out_file_name.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
		if (!this->out_file.is_open())
			throw std::runtime_error("Can't open out file: '" + out_file_name + "'");

		this->out_zip.push(boost::iostreams::gzip_compressor());
		this->out_zip.push(this->out_file);

		if (save_reads_names)
		{
			this->out_reads_file.open(this->reads_file_name.c_str(), std::ios_base::out);
			if (this->out_reads_file.fail())
				throw std::runtime_error("Can't open out reads file: '" + this->reads_file_name + "'");
		}

		this->out_reads_zip.push(boost::iostreams::gzip_compressor());
		this->out_reads_zip.push(this->out_reads_file);
	}

	std::string FilesProcessor::get_out_filename() const
	{
		std::stringstream ss;
		ss << this->base_name << "." << this->out_file_index << ".fastq.gz";
		return ss.str();
	}

	void FilesProcessor::increase_out_file()
	{
		this->out_file_index++;
		std::string out_file_name = this->get_out_filename();

		this->out_zip.pop();
		this->out_file.close();
		this->out_file.open(out_file_name.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
		if (!this->out_file.is_open())
			throw std::runtime_error("Can't open file: '" + out_file_name + "'");

		this->out_zip.push(this->out_file);
	}

	bool FilesProcessor::write(const std::string &text, size_t max_reads)
	{
		bool result = false;
		if (this->current_file_reads_written > max_reads)
		{
			this->increase_out_file();

			this->current_file_reads_written = 0;
			result = true;
		}

		bool out_res = (this->out_zip << text).good();
		if (!out_res)
			throw std::runtime_error("Can't write to file! Text:\n'" + text + "'");

		this->current_file_reads_written++;

		return result;
	}

	void FilesProcessor::write_read_params(const std::string &id_prefix, const Tools::ReadParameters &read_params)
	{
		this->out_reads_zip << read_params.encoded_params(id_prefix) << "\n";
	}

	FilesProcessor::FastQRecord FilesProcessor::get_fastq_record(size_t index)
	{
		FastQRecord record;
		auto &provider = *this->line_providers.at(index);
		if (!provider.get_line(record.id))
			return record;

		if (record.id.at(0) != '@')
			throw std::runtime_error("File '" + this->filenames[index] + "', read '" + record.id + "': fastq malformed!");

		if (!provider.get_line(record.sequence))
			throw std::runtime_error("File '" + this->filenames[index] + "', read '" + record.id + "': fastq ended prematurely!");

		if (!provider.get_line(record.description))
			throw std::runtime_error("File '" + this->filenames[index] + "', read '" + record.id + "': fastq ended prematurely!");

		if (!provider.get_line(record.quality))
			throw std::runtime_error("File '" + this->filenames[index] + "', read '" + record.id + "': fastq ended prematurely!");

		if (record.sequence.length() != record.quality.length())
			throw std::runtime_error("File '" + this->filenames[index] + "', read '" + record.id + "': different lengths of the sequence and the quality string!");

		return record;
	}

	const std::string &FilesProcessor::filename(size_t index) const
	{
		return this->filenames.at(index);
	}

	FilesProcessor::FastQRecord::FastQRecord(const std::string &id, const std::string &sequence,
											 const std::string &description, const std::string &quality)
		: id(id)
		, sequence(sequence)
		, description(description)
		, quality(quality)
	{}
}
