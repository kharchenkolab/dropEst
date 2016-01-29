#include "ResultPrinter.h"

#include "Tools/log_defs.h"

#include <fstream>

void ResultPrinter::print_text_table(const std::string &output_name, const CountMatrix &count_matrix)
{
	L_TRACE << "Writing output matrix to " << output_name << " ";

	std::ofstream output_file(output_name.c_str(), std::ios_base::out);
	// header
	output_file << "gene";
	for (auto const &cell_name: count_matrix.cell_names)
	{
		output_file << '\t' << cell_name;
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

void ResultPrinter::print_binary(const std::string &bin_output_name, const IndropResult &results)
{
	L_TRACE << "writing binary results to " << bin_output_name;

	std::ofstream bin_out_file(bin_output_name.c_str(), std::ios_base::out | std::ios_base::binary);
	boost::archive::binary_oarchive oa(bin_out_file);
	oa << results;
	bin_out_file.close();

	L_TRACE << "All done";
}

void ResultPrinter::print_fields(const std::string &output_suffix, const IndropResult &results)
{
	std::ofstream of("ex_names" + output_suffix);
	for (auto & name : results.exon_count_names)
	{
		of << name << "\n";
	}
	of.close();
	of.open("ex_counts" + output_suffix);
	for (auto & val : results.exon_counts)
	{
		of << val << "\n";
	}
	of.close();
	of.open("merge_n" + output_suffix);
	for (auto & val : results.merge_n)
	{
		of << val << "\n";
	}
	of.close();
	of.open("nonex_names" + output_suffix);
	for (auto & val : results.non_exon_count_names)
	{
		of << val << "\n";
	}
	of.close();
	of.open("nonex_counts" + output_suffix);
	for (auto & val : results.non_exon_counts)
	{
		of << val << "\n";
	}
	of.close();
	of.open("reads_per_umi" + output_suffix);
	for (auto & val : results.reads_per_umi)
	{
		of << val << "\n";
	}
	of.close();
	of.open("umig_covered" + output_suffix);
	for (auto & val : results.umig_covered)
	{
		of << val << "\n";
	}
	of.close();
}
