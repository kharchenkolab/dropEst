#pragma once

//#include <boost/unordered_map.hpp>
#include <string>
#include <unordered_map>

namespace Tools
{
	class ReadsParameters
	{
	public:
		const long total_reads_read;
		const std::string cell_barcode;
		const std::string umi_barcode;

		const bool is_empty;

		template<class Archive>
		void serialize(Archive &ar, const unsigned int /* file_version */)
		{
			ar & this->total_reads_read & this->cell_barcode & this->umi_barcode & this->is_empty;
		}

	public:
		ReadsParameters();
		ReadsParameters(const ReadsParameters &source);
		ReadsParameters(long total_reads_read, const std::string &cell_barcode, const std::string &umi_barcode);
		std::string to_string();
	};



//	typedef boost::unordered_map<std::string, ReadsParameters> reads_params_map_t;
	typedef std::unordered_map<std::string, ReadsParameters> reads_params_map_t;
}