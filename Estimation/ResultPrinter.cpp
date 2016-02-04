#include "ResultPrinter.h"

#include "Tools/log_defs.h"

#include <fstream>
#ifdef R_LIBS
#include <RInside.h>
#endif

void ResultPrinter::print_text_table(const std::string &output_name, const CountMatrix &count_matrix)
{
	L_TRACE << "Writing output matrix to " << output_name << " ";

	std::ofstream output_file(output_name.c_str(), std::ios_base::out);
	// header
	output_file << "gene";
//	for (auto const &cell_name: count_matrix.cell_names)
	for (CountMatrix::s_list_t::const_iterator cell_name_it = count_matrix.cell_names.begin();
		 cell_name_it != count_matrix.cell_names.end(); ++cell_name_it)
	{
		output_file << '\t' << *cell_name_it;
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

	L_TRACE << "Done";
}

void ResultPrinter::print_rds(const std::string &output_name, const IndropResult &results)
{
#ifdef R_LIBS
	L_TRACE << "writing R data to " << output_name;
	RInside R(0, 0);
	R["saved_vec"] = results.get_r_table(output_name);
	R.parseEvalQ("saveRDS(saved_vec, '" + output_name + "')");
	L_TRACE << "Done";
#else
	L_ERR << "Can't print rds without RCpp";
#endif
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
