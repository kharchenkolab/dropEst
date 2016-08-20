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
			L_TRACE << "Writing output matrix to " << output_name << " ";

			std::ofstream output_file(output_name.c_str(), std::ios_base::out);
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

			L_TRACE << "Done";
		}

		void ResultPrinter::print_fields(const std::string &output_suffix, const IndropResult &results)
		{
			throw std::runtime_error("print_fields not implemented");
//	std::ofstream of("ex_names" + output_suffix);
//	for (auto & name : results.exon_chr_count_names)
//	{
//		of << name << "\n";
//	}
//	of.close();
//	of.open("ex_counts" + output_suffix);
//	for (auto & val : results.exon_chr_counts)
//	{
//		of << val << "\n";
//	}
//	of.close();
//	of.open("merge_n" + output_suffix);
//	for (auto & val : results.merge_n)
//	{
//		of << val << "\n";
//	}
//	of.close();
//	of.open("nonex_names" + output_suffix);
//	for (auto & val : results.non_exon_chr_count_names)
//	{
//		of << val << "\n";
//	}
//	of.close();
//	of.open("nonex_counts" + output_suffix);
//	for (auto & val : results.non_exon_chr_counts)
//	{
//		of << val << "\n";
//	}
//	of.close();
//	of.open("reads_per_umi" + output_suffix);
//	for (auto & val : results.reads_per_umi)
//	{
//		of << val << "\n";
//	}
//	of.close();
//	of.open("umig_covered" + output_suffix);
//	for (auto & val : results.umig_covered)
//	{
//		of << val << "\n";
//	}
//	of.close();
		}
	}
}