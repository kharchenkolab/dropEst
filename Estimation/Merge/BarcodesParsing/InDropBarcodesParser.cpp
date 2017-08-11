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

		std::string line;
		Tools::ReverseComplement rc;
		while (std::getline(cb_f, line))
		{
			size_t space_ind = line.find(' ');
			if (space_ind == std::string::npos)
			{
				L_WARN << "WARNING: barcodes line has bad format: '" << line << "'";
				continue;
			}

			barcodes[0].push_back(rc.rc(line.substr(0, space_ind)));
			barcodes[1].push_back(rc.rc(line.substr(space_ind + 1)));
		}

		return barcodes;
	}

	std::vector<std::string> InDropBarcodesParser::split_barcode(const std::string &barcode) const
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
