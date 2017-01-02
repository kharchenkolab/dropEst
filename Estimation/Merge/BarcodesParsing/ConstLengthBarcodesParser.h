#pragma once

#include "BarcodesParser.h"

namespace Estimation
{
	namespace Merge
	{
		namespace BarcodesParsing
		{
			class ConstLengthBarcodesParser : public BarcodesParser
			{
			private:
				std::vector<size_t> _barcode_lengths;
				size_t _barcode_length;

			protected:
				virtual std::vector<std::string> split_barcode(const std::string &barcode) const override;
				virtual barcode_parts_list_t get_barcodes_list(const std::string &barcodes_filename) const override;

			public:
				ConstLengthBarcodesParser(const std::string &barcodes_filename);
				virtual void init() override;
			};
		}
	}
}