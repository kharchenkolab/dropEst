#include <Tools/Logs.h>
#include "ReadMapBamProcessor.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

namespace Estimation
{
	ReadMapBamProcessor::ReadMapBamProcessor(size_t read_prefix_length, const std::string & reads_params_names_str,
	                                         const std::string & gtf_path)
		: BamProcessor(read_prefix_length, gtf_path)
	{
		this->fill_names_map(reads_params_names_str);
	}

	bool ReadMapBamProcessor::get_read_params(const BamTools::BamAlignment &alignment,
	                                          Tools::ReadParameters &read_params) const
	{
		const std::string &read_name = alignment.Name;
		auto iter = this->_reads_params.find(read_name);
		if (iter == this->_reads_params.end())
		{
			L_WARN << "WARNING: can't find read_name in map: " << read_name;
			return false;
		}

		read_params = iter->second;
		if (read_params.is_empty())
		{
			L_WARN << "WARNING: empty parameters for read_name: " << read_name;
			return false;
		}

		return true;
	}

	void ReadMapBamProcessor::fill_names_map(const std::string &reads_params_names_str)
	{
		std::vector<std::string> params_names;
		boost::split(params_names, reads_params_names_str, boost::is_any_of(" \t"));
		size_t total_reads_count = 0;
		L_TRACE << "Start loading reads names";
		for (auto const &name : params_names)
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
				std::string key = row.substr(0, space_index);
				std::string value = row.substr(space_index + 1);
				if (this->_reads_params.find(key) != this->_reads_params.end())
				{
					L_ERR << "Read name is already in map: " << key << ", old value: '" <<
					      this->_reads_params[key].to_monolithic_string() << "', new value: " << value;
					continue;
				}

				try
				{
					this->_reads_params[key] = Tools::ReadParameters(value);
				}
				catch (std::runtime_error err)
				{
					L_ERR << err.what();
					continue;
				}

			}
		}
		L_TRACE << "All reads names were loaded";
	}
}

