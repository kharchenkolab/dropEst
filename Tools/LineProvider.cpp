#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "LineProvider.h"
#include "Logs.h"
#include <zlib.h>
#include <chrono>

namespace Tools
{
	LineProvider::LineProvider(const std::string &filename, size_t buff_size)
		: _buff_size(buff_size)
		, _cur_buff_size(0)
		, _cur_break_index(0)
		, _buff_start(0)
		, _buffer_raw(new char[buff_size + 1])
	{
		this->_in_file = gzopen(filename.c_str(),"rb");
		this->_buffer_raw.get()[buff_size] = 0;
//		if (!this->_in_file) // TODO
//			throw std::runtime_error("Can't open fastq file '" + filename + "'");

//		if (boost::ends_with(filename, ".gz") || boost::ends_with(filename, ".gzip")) // TODO
//		{
//			this->_in_fstream.push(boost::iostreams::gzip_decompressor());
//		}

		gzrewind(this->_in_file);
	}

	bool LineProvider::get_line(std::string &line, bool extension)
	{
		if (this->_buff_start == this->_cur_buff_size)
		{
			if (gzeof(this->_in_file))
				return false;

			this->read_buffer();
		}

		if (this->_cur_break_index == this->_line_breaks.size())
		{
			line = std::string(this->_buffer_raw.get(), this->_buff_start, this->_cur_buff_size - this->_buff_start);
		}
		else
		{
			size_t index = this->_line_breaks[this->_cur_break_index];
			line = std::string(this->_buffer_raw.get(), this->_buff_start, index - this->_buff_start);
		}

		if (line.empty())
		{
			this->_buff_start++;
			this->_cur_break_index++;

			if (extension)
				return true;

			return this->get_line(line);
		}

		if (this->_cur_break_index != this->_line_breaks.size())
		{
			this->_buff_start = this->_line_breaks[this->_cur_break_index] + 1;
			this->_cur_break_index++;
			return true;
		}

		std::string line_extension = "";

		this->read_buffer();
		bool res = this->get_line(line_extension, true);
		line += line_extension;
		return res;
	}

	void LineProvider::read_buffer()
	{
		const char separators[] = "\r\n";

		this->_cur_buff_size = gzread(this->_in_file, this->_buffer_raw.get(), this->_buff_size);
		this->_line_breaks.clear();
		this->_line_breaks.reserve(1000);

		int err;
		auto msg = std::string(gzerror(this->_in_file, &err));

		if (err != 0)
			throw std::runtime_error("File read error: " + msg);

		this->_buff_start = 0;

		char *begin = this->_buffer_raw.get();
		char *c_ptr = begin;

		do {
			char *end = strpbrk(c_ptr, separators);
			if(end)
			{
				this->_line_breaks.push_back(end - begin);
				c_ptr = end + 1;
				continue;
			}

			break;
		} while(*c_ptr);

		this->_cur_break_index = 0;
	}

	LineProvider::~LineProvider()
	{
		gzclose(this->_in_file);
	}
}
