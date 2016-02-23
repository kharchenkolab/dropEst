#include "ReadsParameters.h"

#include <sstream>
#include <string>

namespace Tools
{
	ReadsParameters::ReadsParameters(long total_reads_read, const std::string &cell_barcode, const std::string &umi_barcode)
		: _total_reads_read(total_reads_read)
		, _cell_barcode(cell_barcode)
		, _umi_barcode(umi_barcode)
		, _is_empty(false)
	{
		if (total_reads_read < 0 || cell_barcode.length() == 0 || umi_barcode.length() == 0)
			throw std::runtime_error("Bad reads parameters: " + cell_barcode + umi_barcode);
	}

	ReadsParameters::ReadsParameters(const std::string &encoded_params)
	{
		size_t umi_start_pos = encoded_params.rfind('#');
		if (umi_start_pos == std::string::npos)
			throw std::runtime_error("WARNING: unable to parse out UMI in: " + encoded_params);

		size_t cell_barcode_start_pos = encoded_params.rfind('!', umi_start_pos);
		if (cell_barcode_start_pos == std::string::npos)
			throw std::runtime_error("WARNING: unable to parse out cell tag in: " + encoded_params);

		this->_total_reads_read = atol(encoded_params.substr(0, umi_start_pos).c_str());
		this->_cell_barcode = encoded_params.substr(cell_barcode_start_pos + 1, umi_start_pos - cell_barcode_start_pos - 1);
		this->_umi_barcode = encoded_params.substr(umi_start_pos + 1);
	}

	ReadsParameters::ReadsParameters()
		: _total_reads_read(0)
		, _cell_barcode("")
		, _umi_barcode("")
		, _is_empty(true)
	{}

	ReadsParameters::ReadsParameters(const ReadsParameters &source)
		: _total_reads_read(source._total_reads_read)
		, _cell_barcode(source._cell_barcode)
		, _umi_barcode(source._umi_barcode)
		, _is_empty(source._is_empty)
	{}

	long ReadsParameters::total_reads_read()
	{
		return this->_total_reads_read;
	}

	std::string ReadsParameters::cell_barcode()
	{
		return this->_cell_barcode;
	}

	std::string ReadsParameters::umi_barcode()
	{
		return this->_umi_barcode;
	}

	bool ReadsParameters::is_empty()
	{
		return this->_is_empty;
	}

	std::string ReadsParameters::to_string()
	{
		std::ostringstream text;
		text << '@' << this->_total_reads_read << '!' << this->_cell_barcode << '#' << this->_umi_barcode;
		return text.str();
	}
}