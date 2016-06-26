#include "ReadParameters.h"

#include <sstream>
#include <iostream>
#include <string>

namespace Tools
{
	ReadParameters ReadParameters::from_string(const std::string &str)
	{
		std::istringstream inp(str);

		long read_number;
		std::string read_name;
		std::string cell_barcode;
		std::string umi_barcode;

		inp >> read_number >> read_name >> cell_barcode >> umi_barcode;

		if (inp.fail())
			throw std::runtime_error("Can't parse read string: '" + str + "'");

		return Tools::ReadParameters(read_number, read_name, cell_barcode, umi_barcode);
	}

	std::string ReadParameters::to_string()
	{
		std::ostringstream out;

		out << this->read_number() << " " << this->read_name() << " " <<
				this->cell_barcode() << " " << this->umi_barcode();

		return out.str();
	}

	ReadParameters::ReadParameters(long read_number, const std::string &read_name, const std::string &cell_barcode,
								   const std::string &umi_barcode)
		: _read_number(read_number)
		, _read_name(read_name)
		, _cell_barcode(cell_barcode)
		, _umi_barcode(umi_barcode)
		, _is_empty(false)
	{
		if (read_number < 0 || cell_barcode.length() == 0 || umi_barcode.length() == 0)
			throw std::runtime_error("Bad reads parameters: '" + std::to_string(read_number)  + "' '" + cell_barcode +
											 "' '" + umi_barcode + "'");
	}

	ReadParameters::ReadParameters(const std::string &monolithic_params_string)
		: _is_empty(false)
	{
		size_t start_pos = monolithic_params_string[0] == '@' ? 1 : 0;
		size_t umi_start_pos = monolithic_params_string.rfind('#');
		if (umi_start_pos == std::string::npos)
			throw std::runtime_error("ERROR: unable to parse out UMI in: " + monolithic_params_string);

		size_t cell_barcode_start_pos = monolithic_params_string.rfind('!', umi_start_pos);
		if (cell_barcode_start_pos == std::string::npos)
			throw std::runtime_error("ERROR: unable to parse out cell tag in: " + monolithic_params_string);

		this->_read_number = atol(monolithic_params_string.substr(start_pos, umi_start_pos - start_pos).c_str());
		this->_read_name = monolithic_params_string;
		this->_cell_barcode = monolithic_params_string.substr(cell_barcode_start_pos + 1, umi_start_pos - cell_barcode_start_pos - 1);
		this->_umi_barcode = monolithic_params_string.substr(umi_start_pos + 1);
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

	std::string ReadParameters::to_monolithic_string() const
	{
		std::ostringstream text;
		text << '@' << this->_read_number << '!' << this->_cell_barcode << '#' << this->_umi_barcode;
		return text.str();
	}
}
