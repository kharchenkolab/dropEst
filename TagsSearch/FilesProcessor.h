#pragma once

#include <fstream>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Tools/ReadParameters.h"

namespace TagsSearch
{
	class FilesProcessor
	{
	public:
		struct FastQRecord
		{
			FastQRecord(const std::string &id="", const std::string &sequence="", const std::string &description="",
						const std::string &quality="");

			std::string id;
			std::string sequence;
			std::string description;
			std::string quality;
		};
	private:
		const std::vector<std::string> filenames;
		std::vector<std::shared_ptr<std::ifstream>> in_files;
		std::vector<std::shared_ptr<boost::iostreams::filtering_istream>> in_fstreams;

		std::ofstream out_file;
		std::ofstream out_reads_file;
		boost::iostreams::filtering_ostream out_reads_zip;
		boost::iostreams::filtering_ostream out_zip;

		const std::string base_name;
		const std::string reads_file_name;

		size_t current_file_reads_written;
		size_t out_file_index;

	private:
		std::string get_out_filename() const;
		void increase_out_file();

	public:
		FilesProcessor(const std::vector<std::string> &filenames, const std::string &base_name,
					   bool save_reads_names = false);

		const std::string& filename(size_t index) const;
		FastQRecord get_fastq_record(size_t index);

		bool write(const std::string &text, size_t max_reads);
		void write_read_params(const std::string &id_prefix, const Tools::ReadParameters &read_params);
	};
}
