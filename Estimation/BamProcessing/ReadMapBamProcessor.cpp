#include <Tools/Logs.h>
#include "ReadMapBamProcessor.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

namespace Estimation
{
namespace BamProcessing
{
	ReadMapBamProcessor::ReadMapBamProcessor(size_t read_prefix_length, const std::string & read_param_names,
											 const std::string & gtf_path)
			: BamProcessor(read_prefix_length, gtf_path)
			, _read_param_names(read_param_names)
	{}

	bool ReadMapBamProcessor::get_read_params(const BamTools::BamAlignment &alignment,
											  Tools::ReadParameters &read_params) const
	{
		const std::string &read_name = alignment.Name;
		Tools::reads_params_map_t::const_iterator iter;
		bool read_not_found;
		#pragma omp critical(READ_PARAMS)
		{
			iter = this->_reads_params.find(read_name);
			read_not_found = (iter == this->_reads_params.end());

			if (!read_not_found)
			{
				read_params = iter->second;
				this->_reads_params.erase(iter);
			}
		}

		if (read_not_found)
		{
			L_WARN << "WARNING: can't find read_name in map: " << read_name;
			return false;
		}

		if (read_params.is_empty())
		{
			L_WARN << "WARNING: empty parameters for read_name: " << read_name;
			return false;
		}

		return true;
	}

	void ReadMapBamProcessor::init_temporaries_before_parsing(bool save_read_name) const
	{
		std::vector<std::string> params_names;
		boost::split(params_names, this->_read_param_names, boost::is_any_of(" \t"));
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
					this->_reads_params[key] = Tools::ReadParameters(value, save_read_name);
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

	void ReadMapBamProcessor::release_temporaries_after_parsing() const
	{
		this->_reads_params.clear();
	}
}
}

