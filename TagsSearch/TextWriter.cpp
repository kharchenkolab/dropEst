#include "TextWriter.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace TagsSearch
{
	TextWriter::TextWriter(const std::string &out_file_name, const std::string &file_extension, size_t max_file_size)
		: _out_file_name(out_file_name)
		, _out_file_extension(file_extension)
		, _max_file_size(max_file_size)
		, _limited_file_size(_max_file_size != 0)
		, _current_file_reads_written(0)
		, _out_file_index(0)
	{
		this->increase_out_file();
	}

	bool TextWriter::write(const std::string &text)
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

	std::string TextWriter::get_out_filename() const
	{
		return this->_out_file_name + "." + std::to_string(this->_out_file_index) + "." + this->_out_file_extension;
	}

	std::string TextWriter::gzip(const std::string &text)
	{
		std::string res;
		boost::iostreams::filtering_ostream ofs;
		ofs.push(boost::iostreams::gzip_compressor());
		ofs.push(std::back_inserter(res));

		ofs << text << std::flush;
		return res;
	}

	void TextWriter::increase_out_file()
	{
		this->_out_file_index++;
		std::string out_file_name = this->get_out_filename();

		this->_out_file.close();
		this->_out_file.open(out_file_name.c_str(), std::ios_base::out | std::ios_base::binary);
		if (!this->_out_file.is_open())
			throw std::runtime_error("Can't open file: '" + out_file_name + "'");
	}
}