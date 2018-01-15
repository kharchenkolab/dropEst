#include <Tools/Logs.h>
#include "ReadMapParamsParser.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

namespace Estimation
{
namespace BamProcessing
{
	ReadMapParamsParser::ReadMapParamsParser(const std::string &gtf_path, const std::string &read_param_filenames,
	                                         const BamTags &tags, bool gene_in_chromosome_name, int min_quality)
		: ReadParamsParser(gtf_path, tags, gene_in_chromosome_name)
		, _min_quality(min_quality)
	{
		this->init(read_param_filenames);
	}

	bool ReadMapParamsParser::get_read_params(const BamTools::BamAlignment &alignment,
	                                          Tools::ReadParameters &read_params)
	{
		const std::string &read_name = alignment.Name;

		auto iter = this->_read_params.find(read_name);
		bool read_not_found = (iter == this->_read_params.end());

		if (!read_not_found)
		{
			read_params = iter->second.parameters(this->_barcode_indexer, this->_umi_indexer, this->_umi_quality_indexer);
			this->_read_params.erase(iter);
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

	void ReadMapParamsParser::init(const std::string &read_param_filenames)
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
				if (total_reads_count % 10000000 == 0)
				{
					L_TRACE << "Total " << total_reads_count << " read records processed" << std::flush;
				}

				try
				{
					auto params_pair = Tools::ReadParameters::parse_from_string(row, this->_min_quality);
					auto params = ReadParametersEfficient(params_pair.second, this->_barcode_indexer, this->_umi_indexer,
					                                      this->_umi_quality_indexer);
					auto ins_iter = this->_read_params.emplace(params_pair.first, params);

					if (!ins_iter.second)
					{
						L_ERR << "Read name is already in map: " << ins_iter.first->first << ", old value: '"
						      << ins_iter.first->second.parameters(this->_barcode_indexer, this->_umi_indexer,
						                                           this->_umi_quality_indexer).to_string("").substr(1);
						continue;
					}
				}
				catch (std::runtime_error &err)
				{
					L_ERR << err.what();
					continue;
				}

			}
		}

		L_TRACE << "All read parameters were loaded";
	}
}
}

