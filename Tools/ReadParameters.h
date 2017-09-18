#pragma once

#include <boost/unordered_map.hpp>
#include <string>
#include <map>

namespace Tools
{
	class ReadParameters
	{
	private:
		std::string _cell_barcode;
		std::string _umi;
		std::string _cell_barcode_quality;
		std::string _umi_quality;

		bool _is_empty;

	public:
		static const char quality_offset = 33;

	public:
		ReadParameters();
		ReadParameters(const std::string &cell_barcode, const std::string &umi,
		               const std::string &cell_barcode_quality = "", const std::string &umi_quality = "");

		static ReadParameters parse_encoded_id(const std::string &encoded_id);
		static std::pair<std::string, ReadParameters> parse_from_string(const std::string &params_string);

		std::string to_string(const std::string &read_name) const;
		std::string encoded_id(const std::string &id_prefix) const;

		const std::string& cell_barcode() const;
		const std::string& umi() const;
		const std::string& cell_barcode_quality() const;
		const std::string& umi_quality() const;

		bool is_empty() const;
	};

	typedef boost::unordered_map<std::string, ReadParameters> reads_params_map_t;
//	typedef std::map<std::string, ReadParameters> reads_params_map_t;
}
