#include "ResultPrinter.h"

#include "Tools/Logs.h"
#include <Estimation/Results/IndropResults.h>

#include <fstream>

namespace Estimation
{
	namespace Results
	{
		void ResultPrinter::print_text_table(const std::string &output_name, const CountMatrix &count_matrix)
		{
			Tools::trace_time("Writing output matrix to " + output_name);

			std::ofstream output_file(output_name.c_str(), std::ios_base::out);
			if (output_file.fail())
				throw std::runtime_error("Can't open matrix file: '" + output_name + "'");

			// header
			output_file << "gene";
			for (auto const &cell_barcode: count_matrix.cell_names)
			{
				output_file << '\t' << cell_barcode;
			}
			output_file << "\n";

			size_t rows_count = count_matrix.gene_names.size();
			size_t columns_coult = count_matrix.cell_names.size();
			for (size_t row = 0; row < rows_count; row++)
			{
				output_file << count_matrix.gene_names[row];
				for (size_t col = 0; col < columns_coult; col++)
				{
					output_file << '\t' << count_matrix.counts[row * columns_coult + col];
				}
				output_file << "\n";
			}

			output_file.close();

			Tools::trace_time("Done");
		}
	}
}