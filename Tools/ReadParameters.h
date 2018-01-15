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

		bool _pass_quality_threshold;
		bool _is_empty;

	public:
		static const char quality_offset = 33; // TODO: to xml

	private:
		bool check_quality(int min_quality) const;

	public:
		ReadParameters();
		ReadParameters(const std::string &cell_barcode, const std::string &umi,
		               const std::string &cell_barcode_quality, const std::string &umi_quality,
		               int min_quality);
		ReadParameters(const std::string &cell_barcode, const std::string &umi,
		               const std::string &cell_barcode_quality, const std::string &umi_quality,
		               bool pass_quality_threshold);

		static ReadParameters parse_encoded_id(const std::string &encoded_id);
		static std::pair<std::string, ReadParameters> parse_from_string(const std::string &params_string,
		                                                                int min_quality = 0);

		std::string to_string(const std::string &read_name) const;
		std::string encoded_id(const std::string &id_prefix) const;

		const std::string& cell_barcode() const;
		const std::string& umi() const;
		const std::string& cell_barcode_quality() const;
		const std::string& umi_quality() const;

		bool pass_quality_threshold() const;
		bool is_empty() const;
	};

	using reads_params_map_t = boost::unordered_map<std::string, ReadParameters>;
//	using reads_params_map_t = std::map<std::string, ReadParameters>;
}
