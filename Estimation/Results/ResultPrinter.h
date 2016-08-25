#pragma once

#include <string>

namespace Estimation
{
	namespace Results
	{
		class IndropResult;
		class CountMatrix;

		class ResultPrinter
		{
		private:

		public:
			static void print_text_table(const std::string &output_name, const CountMatrix &count_matrix);

			static void print_fields(const std::string &output_suffix, const IndropResult &results);

		};
	}
}