#pragma once

#include "BarcodesParser.h"

namespace TestEstimator
{
	struct testSplitBarcode;
	struct testConstLengthBarcodeParser;
}

namespace Estimation
{
	namespace Merge
	{
		namespace BarcodesParsing
		{
			class ConstLengthBarcodesParser : public BarcodesParser
			{
				friend struct TestEstimator::testSplitBarcode;
				friend struct TestEstimator::testConstLengthBarcodeParser;

			private:
				std::vector<size_t> _barcode_lengths;
				size_t _barcode_length;

			protected:
				std::vector<std::string> split_barcode(const std::string &barcode) const override;
				barcode_parts_list_t get_barcodes_list(const std::string &barcodes_filename) const override;

			public:
				explicit ConstLengthBarcodesParser(const std::string &barcodes_filename);
				void init() override;
			};
		}
	}
}