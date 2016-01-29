#pragma once

#include "IndropResults.h"

#include <string>

class ResultPrinter
{
private:

public:
	static void print_binary(const std::string &bin_output_name, const IndropResult &results);
	static void print_text_table(const std::string &output_name, const CountMatrix &count_matrix);
	static void print_fields(const std::string &output_suffix, const IndropResult &results);

};

