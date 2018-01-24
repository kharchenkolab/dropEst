#include "InDropBarcodesParser.h"

#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>

#include <fstream>
#include <stdexcept>

namespace Estimation
{
namespace Merge
{
namespace BarcodesParsing
{
	BarcodesParser::barcode_parts_list_t InDropBarcodesParser ::get_barcodes_list(const std::string &barcodes_filename) const
	{
		std::ifstream cb_f(barcodes_filename);
		if (cb_f.fail())
			throw std::runtime_error("Can't open barcodes file: '" + barcodes_filename + "'");

		barcode_parts_list_t barcodes(2);
		for (size_t i = 0; i < 2; ++i)
		{
			if (!BarcodesParser::read_line(cb_f, barcodes[i]) || barcodes[i].empty())
				throw std::runtime_error("File with barcodes (" + barcodes_filename + ") has wrong format");
		}

		return barcodes;
	}

	BarcodesParser::barcodes_list_t InDropBarcodesParser::split_barcode(const std::string &barcode) const
	{
		std::vector<std::string> res;
		res.push_back(barcode.substr(0, barcode.length() - this->_barcode2_length));
		res.push_back(barcode.substr(barcode.length() - this->_barcode2_length));

		return res;
	}

	InDropBarcodesParser::InDropBarcodesParser(const std::string &barcodes_filename)
		: BarcodesParser(barcodes_filename)
	{}

	void InDropBarcodesParser::init()
	{
		BarcodesParser::init();
		this->_barcode2_length = this->barcode(1, 0).length();
	}
}
}
}
