#include "ReadsParameters.h"

#include <sstream>
#include <string>

namespace Tools
{
	ReadsParameters::ReadsParameters(long total_reads_read, const std::string &cell_barcode, const std::string &umi_barcode)
		: total_reads_read(total_reads_read)
		, cell_barcode(cell_barcode)
		, umi_barcode(umi_barcode)
		, is_empty(false)
	{
		if (total_reads_read < 0 || cell_barcode.length() == 0 || umi_barcode.length() == 0)
			throw std::runtime_error("Bad reads parameters: " + cell_barcode + umi_barcode);
	}

	std::string ReadsParameters::to_string()
	{
		std::ostringstream text;
		text << '@' << this->total_reads_read << '!' << this->cell_barcode << '#' << this->umi_barcode;
		return text.str();
	}

	ReadsParameters::ReadsParameters()
		: total_reads_read(0)
		, cell_barcode("")
		, umi_barcode("")
		, is_empty(true)
	{}

	ReadsParameters::ReadsParameters(const ReadsParameters &source)
		: total_reads_read(source.total_reads_read)
		, cell_barcode(source.cell_barcode)
		, umi_barcode(source.umi_barcode)
		, is_empty(source.is_empty)
	{}
}