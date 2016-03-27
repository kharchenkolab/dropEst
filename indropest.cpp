#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>
#include <ctime>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Estimation/Estimator.h"
#include "Estimation/Results/ResultPrinter.h"
#include "Estimation/Results/IndropResults.h"
#include "Estimation/Results/BadCellsStats.h"
#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"

using namespace std;
using namespace Estimation;

struct Params
{
	bool bam_output = false;
	bool cant_parse = false;
	bool verbose = false;
	bool text_output = false;
	bool merge_tags = false;
	bool not_filtered = false;
	string output_name = "";
	string config_file_name = "";
	string log_prefix = "";
	string reads_params_file = "";
	string gtf_filename = "";
};

static void usage()
{
	cerr << "\tindropest: estimate molecular counts per cell" << endl;
	cerr << "SYNOPSIS\n";
	cerr <<
	"\tindropest [-t, --text-output] [-m|--merge-cell-tags] [-v|--verbose] [-n | --not-filtered] [-g | --gtf filename]"
	"[-l, --log-prefix logs_name] [-r, --reads-params filename] -c config.xml file1.bam [file2.bam ...] [-b | --bam]"
	<< endl;
	cerr << "OPTIONS:\n";
	cerr << "\t-b, --bam-out: print corrected bam files" << endl;
	cerr << "\t-c, --config filename: xml file with estimation parameters" << endl;
	cerr << "\t-g, --gtf filename: gtf file with genes annotations" << endl;
	cerr << "\t-l, --log-prefix : logs prefix" << endl;
	cerr << "\t-m, --merge-cell-tags : merge linked cell tags" << endl;
	cerr << "\t-n, --not-filtered : print data for all cells" << endl;
	cerr << "\t-o, --output-file filename : output file name" << endl;
	cerr << "\t-r, --reads-params filename: file with serialized params from tags search step" << endl;
	cerr << "\t-t, --text-output : write out text matrix" << endl;
	cerr << "\t-v, --verbose : verbose mode" << endl;
}

static Params parse_cmd_params(int argc, char **argv)
{
	Params params;

	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"bam-output",     	no_argument, 	   0, 'b'},
			{"config",     		required_argument, 0, 'c'},
			{"gtf",     		required_argument, 0, 'g'},
			{"log-prefix",		required_argument, 0, 'l'},
			{"merge-cell-tags", no_argument,       0, 'm'},
			{"not-filtered",	no_argument, 	   0, 'n'},
			{"output-file",     required_argument, 0, 'o'},
			{"reads-params",     required_argument, 0, 'r'},
			{"text-output",     no_argument,       0, 't'},
			{"verbose",         no_argument,       0, 'v'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "bvtmno:c:l:r:g:", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'b':
				params.bam_output = true;
				break;
			case 'v' :
				params.verbose = true;
				break;
			case 'm' :
				params.merge_tags = true;
				break;
			case 't' :
				params.text_output = true;
				break;
			case 'n' :
				params.not_filtered = true;
				break;
			case 'o' :
				params.output_name = string(optarg);
				break;
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			case 'r' :
				params.reads_params_file = string(optarg);
				break;
			case 'l' :
				params.log_prefix = string(optarg);
				break;
			case 'g' :
				params.gtf_filename = string(optarg);
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

	if (params.log_prefix.length() != 0)
	{
		params.log_prefix += "_";
	}
	Tools::init_log(params.verbose, false, params.log_prefix + "est_main.log", params.log_prefix + "est_debug.log");

	vector<string> files;
	while (optind < argc)
	{
		files.push_back(string(argv[optind++]));
	}

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	time_t ctt = time(0);
	try
	{
		L_TRACE << "Run: " << asctime(localtime(&ctt));
		Estimator estimator(pt.get_child("config.Estimation"));
		CellsDataContainer container = estimator.get_cells_container(files, params.merge_tags, params.bam_output,
		                                                             params.reads_params_file, params.gtf_filename);
		
		if (false)
		{
			Results::IndropResult results = estimator.get_results(container, params.not_filtered);
		
			ctt = time(0);
			L_TRACE << "Done: " << asctime(localtime(&ctt));

			if (params.text_output)
			{
				Results::ResultPrinter::print_text_table(params.output_name, results.cm);
				params.output_name += ".bin";
			}

#ifdef R_LIBS
			results.save_rds(params.output_name);
#else
			Results::ResultPrinter::print_binary(params.output_name, results);
#endif
		}
#ifdef R_LIBS
		L_TRACE << "Get bad cells results";
		Results::BadCellsStats bad_cells_results(container.stats());
		L_TRACE << "Done";
		bad_cells_results.save_rds(params.output_name + ".bc.rds");
#endif
	}
	catch (std::runtime_error err)
	{
		L_ERR << err.what();
		return 1;
	}
	ctt = time(0);
	L_TRACE << "All done: " << asctime(localtime(&ctt));
}
