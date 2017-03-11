#pragma once

#include <boost/unordered_map.hpp>
#include <string>
#include <map>

namespace Tools
{
	class ReadParameters
	{
	private:
		[[deprecated("Read names will be removed")]]
		std::string _read_name;

		std::string _cell_barcode;
		std::string _umi;
		std::string _cell_barcode_quality;
		std::string _umi_quality;

		bool _is_empty;

	public:
		ReadParameters();
		ReadParameters(const std::string &monolithic_params_string, bool save_read_name=true);
		ReadParameters(const ReadParameters &source);
		ReadParameters(const std::string &cell_barcode, const std::string &umi, const std::string &cell_barcode_quality,
		               const std::string &umi_quality);

		[[deprecated("Read names will be removed")]]
		ReadParameters(const std::string &read_name, const std::string &cell_barcode, const std::string &umi);
		[[deprecated("Read names will be removed. Use encoded_id instead.")]]
		std::string to_monolithic_string(const std::string &fake_name="") const;
		[[deprecated("Read names will be removed")]]
		std::string read_name_safe() const;
		[[deprecated("Read names will be removed")]]
		const std::string& read_name() const;

		std::string encoded_id(const std::string &id_prefix) const;
		std::string encoded_params(const std::string &id_prefix) const;

		const std::string& cell_barcode() const;
		const std::string& umi() const;
		const std::string& cell_barcode_quality() const;
		const std::string& umi_quality() const;

		bool is_empty() const;
	};

	typedef boost::unordered_map<std::string, ReadParameters> reads_params_map_t;
//	typedef std::map<std::string, ReadParameters> reads_params_map_t;
}