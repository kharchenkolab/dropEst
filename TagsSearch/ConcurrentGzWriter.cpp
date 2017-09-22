#include "ConcurrentGzWriter.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace TagsSearch
{
	const size_t ConcurrentGzWriter::max_cache_size = 50000;

	ConcurrentGzWriter::ConcurrentGzWriter(const std::string &out_file_name, const std::string &file_extension,
	                                       size_t max_file_size)
		: _out_file_name(out_file_name)
		, _out_file_extension(file_extension)
		, _max_file_size(max_file_size)
		, _limited_file_size(_max_file_size != 0)
		, _write_in_progress(false)
		, _lines(ConcurrentGzWriter::max_cache_size)
		, _gzipped(ConcurrentGzWriter::max_cache_size)
		, _current_file_reads_written(0)
		, _out_file_index(0)
	{
		this->increase_out_file();
	}

	bool ConcurrentGzWriter::write(const std::string &text)
	{
		bool result = false;
		if (this->_limited_file_size && this->_current_file_reads_written > this->_max_file_size)
		{
			this->increase_out_file();

			this->_current_file_reads_written = 0;
			result = true;
		}

		if (!(this->_out_file << text))
			throw std::runtime_error("Can't write to file: " + this->get_out_filename());

		this->_current_file_reads_written++;

		return result;
	}

	std::string ConcurrentGzWriter::get_out_filename() const
	{
		std::string name = this->_out_file_name;

		if (this->_limited_file_size)
		{
			name += "." + std::to_string(this->_out_file_index);
		}
		return name + "." + this->_out_file_extension;
	}

	std::string ConcurrentGzWriter::gzip(const std::string &text)
	{
		std::string res;
		boost::iostreams::filtering_ostream ofs;
		ofs.push(boost::iostreams::gzip_compressor());
		ofs.push(std::back_inserter(res));

		ofs << text << std::flush;
		return res;
	}

	void ConcurrentGzWriter::increase_out_file()
	{
		this->_out_file_index++;
		std::string out_file_name = this->get_out_filename();

		this->_out_file.close();
		this->_out_file.open(out_file_name.c_str(), std::ios_base::out | std::ios_base::binary);
		if (!this->_out_file.is_open())
			throw std::runtime_error("Can't open file: '" + out_file_name + "'");
	}

	bool ConcurrentGzWriter::full() const
	{
		return this->_lines.full();
	}

	bool ConcurrentGzWriter::empty() const
	{
		return this->_lines.empty() && this->_gzipped.empty();
	}

	void ConcurrentGzWriter::flush_gzip(bool unlimited_size)
	{
		while (unlimited_size || !this->_gzipped.full())
		{
			std::string line;
			if (!this->_lines.pop(line))
				break;

			this->_gzipped.push(ConcurrentGzWriter::gzip(line));
		}
	}

	void ConcurrentGzWriter::flush_write()
	{
		if (this->_write_in_progress.exchange(true))
			return;

		while (true)
		{
			std::string to_write, cur_text;
			while (to_write.length() < ConcurrentGzWriter::max_cache_size && this->_gzipped.pop(cur_text))
			{
				to_write += cur_text;
			}
			if (to_write.empty())
				break;

			this->write(to_write);
		}

		this->_write_in_progress = false;
	}

	void ConcurrentGzWriter::enqueue_line(const std::string &line)
	{
		this->_lines.push(line);
	}

	const std::string &ConcurrentGzWriter::base_filename() const
	{
		return this->_out_file_name;
	}
}