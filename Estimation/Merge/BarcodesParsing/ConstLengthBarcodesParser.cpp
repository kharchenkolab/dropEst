#include "ConstLengthBarcodesParser.h"

#include <Tools/Logs.h>

#include <fstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>


namespace Estimation
{
namespace Merge
{
namespace BarcodesParsing
{
	ConstLengthBarcodesParser::ConstLengthBarcodesParser(const std::string &barcodes_filename)
		: BarcodesParser(barcodes_filename)
		, _barcode_length(0)
	{}

	void ConstLengthBarcodesParser::init()
	{
		BarcodesParser::init();
		for (size_t i = 0; i < this->barcode_parts_num(); ++i)
		{
			size_t cur_len = this->barcode(i, 0).length();
			this->_barcode_lengths.push_back(cur_len);
			this->_barcode_length += cur_len;
		}
	}

	std::vector<std::string> ConstLengthBarcodesParser::split_barcode(const std::string &barcode) const
	{
		if (barcode.length() != this->_barcode_length)
			throw std::runtime_error("Barcode '" + barcode + "' has wrong length (" + std::to_string(this->_barcode_length) + " expected)");

		std::vector<std::string> res;
		size_t prev_ind = 0;
		for (size_t length : this->_barcode_lengths)
		{
			res.push_back(barcode.substr(prev_ind, length));
			prev_ind += length;
		}

		return res;
	}

	BarcodesParser::barcode_parts_list_t ConstLengthBarcodesParser::get_barcodes_list(const std::string &barcodes_filename) const
	{
		std::ifstream cb_f(barcodes_filename);
		if (cb_f.fail())
			throw std::runtime_error("Can't open barcodes file: '" + barcodes_filename + "'");

		barcode_parts_list_t barcodes;

		barcodes_list_t barcode_part;
		while (BarcodesParser::read_line(cb_f, barcode_part, true))
		{
			if (barcode_part.empty())
				throw std::runtime_error("File with barcodes (" + barcodes_filename + ") has wrong format");

			barcodes.emplace_back(barcode_part);
			barcode_part.clear();
		}

		return barcodes;
	}
}
}
}