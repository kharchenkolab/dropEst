#pragma once

//#include <boost/unordered_map.hpp>
#include <string>
#include <unordered_map>
#include <boost/serialization/access.hpp>

namespace Tools
{
	class ReadsParameters
	{
	private:
		long _total_reads_read;
		std::string _cell_barcode;
		std::string _umi_barcode;

		bool _is_empty;

	private:
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive &ar, const unsigned int /* file_version */)
		{
			ar & this->_total_reads_read & this->_cell_barcode & this->_umi_barcode & this->_is_empty;
		}

	public:
		ReadsParameters();
		ReadsParameters(const std::string &params_string);
		ReadsParameters(const ReadsParameters &source);
		ReadsParameters(long total_reads_read, const std::string &cell_barcode, const std::string &umi_barcode);
		std::string to_string();

		long total_reads_read();
		std::string cell_barcode();
		std::string umi_barcode();

		bool is_empty();
	};

//	typedef boost::unordered_map<std::string, ReadsParameters> reads_params_map_t;
	typedef std::unordered_map<std::string, ReadsParameters> reads_params_map_t;
}