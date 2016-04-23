#pragma once

#include <boost/unordered_map.hpp>
#include <string>
#include <map>
#include <boost/serialization/access.hpp>

namespace Tools
{
	class ReadParameters
	{
	private:
		long _read_number;
		std::string _read_name;
		std::string _cell_barcode;
		std::string _umi_barcode;

		bool _is_empty;

	private:
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive &ar, const unsigned int /* file_version */)
		{
			ar & this->_read_number & this->_cell_barcode & this->_umi_barcode & this->_is_empty;
		}

	public:
		static ReadParameters from_string(const std::string &str);

		ReadParameters();
		ReadParameters(const std::string &monolithic_params_string);
		ReadParameters(const ReadParameters &source);
		ReadParameters(long read_number, const std::string &read_name, const std::string &cell_barcode,
					   const std::string &umi_barcode);
		std::string to_monolithic_string() const;
		std::string to_string();

		long read_number() const;
		const std::string& read_name() const;
		const std::string& cell_barcode() const;
		const std::string& umi_barcode() const;

		bool is_empty() const;
	};

	typedef boost::unordered_map<std::string, ReadParameters> reads_params_map_t;
//	typedef std::map<std::string, ReadParameters> reads_params_map_t;
}