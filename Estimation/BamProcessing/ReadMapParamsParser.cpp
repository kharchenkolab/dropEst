#include <Tools/Logs.h>
#include "ReadMapParamsParser.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

namespace Estimation
{
namespace BamProcessing
{
	ReadMapParamsParser::ReadMapParamsParser(const std::string &gtf_path, bool save_read_names,
											 const std::string &read_param_filenames, const BamTags &tags,
											 bool gene_in_chromosome_name)
		: ReadParamsParser(gtf_path, tags, gene_in_chromosome_name)
	{
		this->init(read_param_filenames, save_read_names);
	}

	bool ReadMapParamsParser::get_read_params(const BamTools::BamAlignment &alignment,
											  Tools::ReadParameters &read_params)
	{
		const std::string &read_name = alignment.Name;

		auto iter = this->_reads_params.find(read_name);
		bool read_not_found = (iter == this->_reads_params.end());

		if (!read_not_found)
		{
			read_params = iter->second;
			this->_reads_params.erase(iter);
		}

		if (read_not_found)
		{
			L_WARN << "WARNING: can't find read name: " << read_name;
			return false;
		}

		if (read_params.is_empty())
		{
			L_WARN << "WARNING: empty parameters for read name: " << read_name;
			return false;
		}

		return true;
	}

	void ReadMapParamsParser::init(const std::string &read_param_filenames, bool save_read_name)
	{
		std::vector<std::string> param_filenames;
		boost::split(param_filenames, read_param_filenames, boost::is_any_of(" \t"));
		size_t total_reads_count = 0;
		L_TRACE << "Start loading read parameters";
		for (auto const &name : param_filenames)
		{
			if (name.empty())
				continue;
			L_TRACE << "Start reading file: " << name;
			std::ifstream ifs(name);
			boost::iostreams::filtering_istream gz_fs;
			gz_fs.push(boost::iostreams::gzip_decompressor());
			gz_fs.push(ifs);
			std::string row;
			while (std::getline(gz_fs, row))
			{
				if (row.empty())
					continue;

				total_reads_count++;
				if (total_reads_count % 1000000 == 0)
				{
					L_TRACE << "Total " << total_reads_count << " reads record processed";
				}

				size_t space_index = row.find_first_of(" \t");
				if (space_index == std::string::npos)
				{
					L_ERR << "Can't parse row with reads params: '" << row << "'";
					continue;
				}

				try
				{
					auto read_params_info = Tools::ReadParameters::parse_from_string(row);
					if (this->_reads_params.find(read_params_info.first) != this->_reads_params.end())
					{
						L_ERR << "Read name is already in map: " << read_params_info.first << ", old value: '"
						      << this->_reads_params[read_params_info.first].to_string("").substr(1)
						      << "', new value: " << read_params_info.second.to_string("").substr(1);
						continue;
					}

					this->_reads_params.insert(read_params_info);
				}
				catch (std::runtime_error err)
				{
					L_ERR << err.what();
					continue;
				}

			}
		}

		L_TRACE << "All read parameters were loaded";
	}

	ReadMapParamsParser::~ReadMapParamsParser()
	{
		this->_reads_params.clear();
	}
}
}

