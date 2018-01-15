#include "ReadParameters.h"

namespace Tools
{
	ReadParameters::ReadParameters(const std::string &cell_barcode, const std::string &umi,
	                               const std::string &cell_barcode_quality, const std::string &umi_quality,
	                               int min_quality)
		: _cell_barcode(cell_barcode)
		, _umi(umi)
		, _cell_barcode_quality(cell_barcode_quality)
		, _umi_quality(umi_quality)
		, _pass_quality_threshold(this->check_quality(min_quality))
		, _is_empty(false)
	{
		if (cell_barcode.length() == 0 || umi.length() == 0)
			throw std::runtime_error("Wrong read parameters: '" + cell_barcode + "' '" + umi + "'");
	}

	ReadParameters::ReadParameters(const std::string &cell_barcode, const std::string &umi,
	                               const std::string &cell_barcode_quality, const std::string &umi_quality,
	                               bool pass_quality_threshold)
		: _cell_barcode(cell_barcode)
		, _umi(umi)
		, _cell_barcode_quality(cell_barcode_quality)
		, _umi_quality(umi_quality)
		, _pass_quality_threshold(pass_quality_threshold)
		, _is_empty(false)
	{
		if (cell_barcode.length() == 0 || umi.length() == 0)
			throw std::runtime_error("Wrong read parameters: '" + cell_barcode + "' '" + umi + "'");
	}

	ReadParameters::ReadParameters()
		: _cell_barcode("")
		, _umi("")
		, _cell_barcode_quality("")
		, _umi_quality("")
		, _pass_quality_threshold(false)
		, _is_empty(true)
	{}

	ReadParameters ReadParameters::parse_encoded_id(const std::string &encoded_id)
	{
		size_t umi_start_pos = encoded_id.rfind('#');
		if (umi_start_pos == std::string::npos)
			throw std::runtime_error("ERROR: unable to parse out UMI in: " + encoded_id);

		size_t cell_barcode_start_pos = encoded_id.rfind('!', umi_start_pos);
		if (cell_barcode_start_pos == std::string::npos)
			throw std::runtime_error("ERROR: unable to parse out cell barcode in: " + encoded_id);

		auto cell_barcode = encoded_id.substr(cell_barcode_start_pos + 1, umi_start_pos - cell_barcode_start_pos - 1);
		auto umi = encoded_id.substr(umi_start_pos + 1);

		return ReadParameters(cell_barcode, umi, "", "", 0);
	}

	std::pair<std::string, ReadParameters> ReadParameters::parse_from_string(const std::string &params_string,
	                                                                         int min_quality)
	{
		std::vector<std::string> parsed;
		std::string::size_type start_pos = 0, end_pos = 0;
		for (int i = 0; i < 4; ++i)
		{
			end_pos = params_string.find(' ', start_pos);
			if (end_pos == std::string::npos)
				throw std::runtime_error("ERROR: can't parse read parameters from string: '" + params_string + "'");

			parsed.push_back(params_string.substr(start_pos, end_pos - start_pos));
			start_pos = end_pos + 1;
		}

		parsed.push_back(params_string.substr(start_pos));

		return std::make_pair(parsed[0], ReadParameters(parsed[1], parsed[2], parsed[3], parsed[4], min_quality));
	}

	const std::string& ReadParameters::cell_barcode() const
	{
		return this->_cell_barcode;
	}

	const std::string& ReadParameters::umi() const
	{
		return this->_umi;
	}

	const std::string& ReadParameters::cell_barcode_quality() const
	{
		return this->_cell_barcode_quality;
	}

	const std::string& ReadParameters::umi_quality() const
	{
		return this->_umi_quality;
	}

	bool ReadParameters::is_empty() const
	{
		return this->_is_empty;
	}

	std::string ReadParameters::to_string(const std::string &read_name) const
	{
		return read_name + ' ' + this->_cell_barcode + ' ' + this->_umi + ' ' + this->_cell_barcode_quality + ' ' +
				this->_umi_quality;
	}

	std::string ReadParameters::encoded_id(const std::string &id_prefix) const
	{
		return id_prefix + '!' + this->_cell_barcode + '#' + this->_umi;
	}

	bool ReadParameters::check_quality(int min_quality) const
	{
		if (min_quality <= 0)
			return true;

		for (char qual : this->_cell_barcode_quality)
		{
			if (qual < min_quality + Tools::ReadParameters::quality_offset)
				return false;
		}

		for (char qual : this->_umi_quality)
		{
			if (qual < min_quality + Tools::ReadParameters::quality_offset)
				return false;
		}

		return true;
	}

	bool ReadParameters::pass_quality_threshold() const
	{
		return this->_pass_quality_threshold;
	}
}
