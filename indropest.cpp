#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Estimation/Estimator.h"
#include "Estimation/ResultPrinter.h"
#include "Tools/log_init.h"
#include "Tools/log_defs.h"

using namespace std;

struct Params
{
	bool cant_parse;
	bool verbose;
	bool text_output;
	bool merge_tags;
	string output_name;
	string config_file_name;

	Params() : cant_parse(false), verbose(false), text_output(false), merge_tags(false)
			, output_name(""), config_file_name("")
	{}
};

static void usage()
{
	cerr << "\tindropest: estimate molecular counts per cell" << endl;
	cerr << "SYNOPSIS\n";
	cerr <<
	"\tindropest [-g|--min-cells_genes 1000] [-u|--min-umis 10000] [-m|--merge-cell-tags] [-R|--output-r] [-v|--verbose] -c config.xml file1.bam [file2.bam ...]" <<
	endl;
	cerr << "OPTIONS:\n";
	cerr << "\t-o, --output-file filename : output file name" << endl;
	cerr << "\t-t, --text-output : write out text matrix" << endl;
	cerr << "\t-c, --config config.xml : xml file with estimation parameters" << endl;
	cerr << "\t-m, --merge-cell-tags : merge linked cell tags" << endl;
	cerr << "\t-R, --otuput-r : write out RData file" << endl;
}

static Params parse_cmd_params(int argc, char **argv)
{
	Params params;

	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"verbose",         no_argument,       0, 'v'},
			{"text-output",     no_argument,       0, 't'},
			{"merge-cell-tags", no_argument,       0, 'm'},
			{"output-file",     required_argument, 0, 'o'},
			{"config",     		required_argument, 0, 'c'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "vtmo:c:", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'v' :
				params.verbose = true;
				break;
			case 'm' :
				params.merge_tags = true;
				break;
			case 't' :
				params.text_output = true;
				break;
			case 'o' :
				params.output_name = string(optarg);
				break;
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			default:
				cerr << "indropest: unknown arguments passed" << endl;
				params.cant_parse = true;
				return params;
		}
	}

	if (optind > argc - 1)
	{
		cerr << "indropset: at least one bam file must be supplied" << endl;
		params.cant_parse = true;
		return params;
	}

	if (params.config_file_name == "")
	{
		cerr << "indropset: config file must be supplied" << endl;
		params.cant_parse = true;
		return params;
	}

	if (params.output_name == "")
	{
		if (params.text_output)
		{
			params.output_name = "cell.counts.txt";
		}
		else
		{
			params.output_name = "cell.counts.bin";
		}
	}

	return params;
}

int main(int argc, char **argv)
{
	Params params = parse_cmd_params(argc, argv);

	if (params.cant_parse)
	{
		usage();
		return 1;
	}

	init_log(params.verbose, true, "main_est.log", "debug_est.log");

	vector<string> files;
	while (optind < argc)
	{
		files.push_back(string(argv[optind++]));
	}

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	try
	{
		Estimator estimator(pt.get_child("config.Estimation"));
		IndropResult results = estimator.get_results(files, params.merge_tags);
		L_TRACE << "Done";

		if (params.text_output)
		{
			ResultPrinter::print_text_table(params.output_name, results.cm);
			params.output_name += ".bin";
		}
#ifdef R_LIBS
		ResultPrinter::print_rds(params.output_name, results);
#else
		ResultPrinter::print_binary(params.output_name, results);
#endif
	}
	catch (std::runtime_error err)
	{
		L_ERR << err.what();
		return 1;
	}
	L_TRACE << "All done";
}
