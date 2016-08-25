#pragma once

#include <boost/unordered_map.hpp>
#include <string>
#include <map>

namespace Tools
{
	class ReadParameters
	{
	private:
		std::string _read_name;
		std::string _cell_barcode;
		std::string _umi_barcode;

		bool _is_empty;

	public:
		ReadParameters();
		ReadParameters(const std::string &monolithic_params_string);
		ReadParameters(const ReadParameters &source);
		ReadParameters(const std::string &read_name, const std::string &cell_barcode, const std::string &umi_barcode);
		std::string to_monolithic_string() const;

		const std::string& read_name() const;
		const std::string& cell_barcode() const;
		const std::string& umi_barcode() const;

		bool is_empty() const;
	};

	typedef boost::unordered_map<std::string, ReadParameters> reads_params_map_t;
//	typedef std::map<std::string, ReadParameters> reads_params_map_t;
}