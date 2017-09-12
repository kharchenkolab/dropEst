#pragma once

#include <fstream>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Tools/ReadParameters.h"

namespace TagsSearch
{
	class FastQReader
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

			std::string to_string() const;
		};
	private:
		const std::string _filename;
		std::ifstream _in_file;
		boost::iostreams::filtering_istream _in_fstream;


	public:
		FastQReader(const std::string &filename);

		const std::string& filename() const;
		FastQRecord get_next_record();
	};
}
