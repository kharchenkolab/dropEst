#include "ReadParameters.h"

#include <sstream>
#include <string>

namespace Tools
{
	ReadParameters::ReadParameters(long total_reads_read, const std::string &read_name, const std::string &cell_barcode,
								   const std::string &umi_barcode)
		: _read_number(total_reads_read)
		, _read_name(read_name)
		, _cell_barcode(cell_barcode)
		, _umi_barcode(umi_barcode)
		, _is_empty(false)
	{
		if (total_reads_read < 0 || cell_barcode.length() == 0 || umi_barcode.length() == 0)
			throw std::runtime_error("Bad reads parameters: " + cell_barcode + umi_barcode);
	}

	ReadParameters::ReadParameters(const std::string &encoded_params)
		: _is_empty(false)
	{
		size_t start_pos = encoded_params[0] == '@' ? 1 : 0;
		size_t umi_start_pos = encoded_params.rfind('#');
		if (umi_start_pos == std::string::npos)
			throw std::runtime_error("WARNING: unable to parse out UMI in: " + encoded_params);

		size_t cell_barcode_start_pos = encoded_params.rfind('!', umi_start_pos);
		if (cell_barcode_start_pos == std::string::npos)
			throw std::runtime_error("WARNING: unable to parse out cell tag in: " + encoded_params);

		this->_read_number = atol(encoded_params.substr(start_pos, umi_start_pos - start_pos).c_str());
		this->_read_name = encoded_params;
		this->_cell_barcode = encoded_params.substr(cell_barcode_start_pos + 1, umi_start_pos - cell_barcode_start_pos - 1);
		this->_umi_barcode = encoded_params.substr(umi_start_pos + 1);
	}

	ReadParameters::ReadParameters()
		: _read_number(0)
		, _read_name("")
		, _cell_barcode("")
		, _umi_barcode("")
		, _is_empty(true)
	{}

	ReadParameters::ReadParameters(const ReadParameters &source)
		: _read_number(source._read_number)
		, _read_name(source._read_name)
		, _cell_barcode(source._cell_barcode)
		, _umi_barcode(source._umi_barcode)
		, _is_empty(source._is_empty)
	{}

	long ReadParameters::read_number() const
	{
		return this->_read_number;
	}

	const std::string& ReadParameters::read_name() const
	{
		return this->_read_name;
	}

	const std::string& ReadParameters::cell_barcode() const
	{
		return this->_cell_barcode;
	}

	const std::string& ReadParameters::umi_barcode() const
	{
		return this->_umi_barcode;
	}

	bool ReadParameters::is_empty() const
	{
		return this->_is_empty;
	}

	std::string ReadParameters::to_string() const
	{
		std::ostringstream text;
		text << '@' << this->_read_number << '!' << this->_cell_barcode << '#' << this->_umi_barcode;
		return text.str();
	}
}
