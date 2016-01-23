#pragma once

#include <fstream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

class FilesProcessor
{
private:
	std::ifstream r1_file;
	std::ifstream r2_file;
	boost::iostreams::filtering_istream r1_fs;
	boost::iostreams::filtering_istream r2_fs;

	std::ofstream out_file;
	boost::iostreams::filtering_ostream out_zip;

	std::string base_name;
	long max_reads;
	int current_file_reads_written;
	int out_file_index;

private:
	std::string get_out_filename();
	void increase_out_file();

public:
	FilesProcessor(const std::string &r1_filename, const std::string &r2_filename,
				   const std::string &base_name, long max_reads);

	bool get_r1_line(std::string &out);
	bool get_r2_line(std::string &out);

	bool write(const std::string &text);

	void close();
};