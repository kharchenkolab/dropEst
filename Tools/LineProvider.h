#pragma once

#include <fstream>
#include <string>

#include <zlib.h>
#include <queue>

namespace Tools
{
	class LineProvider
	{
	private:
		const size_t _buff_size;
		size_t _cur_buff_size;
		size_t _cur_break_index;
		size_t _buff_start;
		std::shared_ptr<char> _buffer_raw;

		gzFile _in_file;
		std::vector<size_t> _line_breaks;

	public:
		LineProvider(const std::string &filename, size_t buff_size);
		virtual ~LineProvider();
		bool get_line(std::string &line, bool extension=false);

		void read_buffer();
	};
}
