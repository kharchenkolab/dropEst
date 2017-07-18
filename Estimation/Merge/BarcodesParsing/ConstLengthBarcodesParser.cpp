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
		std::vector<size_t> part_lengths; //to check that every line has equal lengths

		std::string line;
		while (std::getline(cb_f, line))
		{
			if (line.empty())
				continue;

			std::vector<std::string> cb_parts;
			boost::trim_if(line, boost::is_any_of(" \t\r\n"));
			boost::split(cb_parts, line, boost::is_any_of(" \t\r\n"), boost::token_compress_on);

			if (!part_lengths.empty() && cb_parts.size() != part_lengths.size())
			{
				L_WARN << "WARNING: barcode line has bad format: '" << line << "' (" << cb_parts.size() << " parts, "
					   << part_lengths.size() << " expected)";
				continue;
			}

			if (part_lengths.empty())
			{
				part_lengths.resize(cb_parts.size(), 0);
				barcodes.resize(cb_parts.size());
			}

			for (size_t part_ind = 0; part_ind < part_lengths.size(); ++part_ind)
			{
				if (part_lengths[part_ind] == 0)
				{
					part_lengths[part_ind] = cb_parts[part_ind].length();
				}

				if (part_lengths[part_ind] != cb_parts[part_ind].length())
				{
					L_WARN << "WARNING: barcode part has wrong length: '" << cb_parts[part_ind] << "' ("
						   << this->_barcode_lengths[part_ind] << " expected)";
					continue;
				}
				barcodes[part_ind].push_back(cb_parts[part_ind]);
			}
		}

		return barcodes;
	}
}
}
}