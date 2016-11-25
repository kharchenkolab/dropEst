#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>
#include <ctime>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <omp.h>
#include <Estimation/BamProcessing/BamController.h>

#include "Estimation/Estimator.h"
#include "Estimation/Results/ResultPrinter.h"
#include "Estimation/Results/IndropResults.h"
#include "Estimation/Results/BadCellsStats.h"
#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"
#include "Tools/UtilFunctions.h"

using namespace std;
using namespace Estimation;

struct Params
{
	bool bam_output = false;
	bool cant_parse = false;
	bool filled_bam = false;
	bool filtered_bam_output = false;
	bool merge_tags = false;
	bool not_filtered = false;
	bool reads_output = false;
	bool text_output = false;
	bool verbose = false;
	string barcodes_filename = "";
	string config_file_name = "";
	string genesets_rds = "";
	string gtf_filename = "";
	string log_prefix = "";
	string output_name = "";
	string reads_params_names_str = "";
	int num_of_threads = 1;
};

static void check_files_existence(const Params &params, const vector<string> &bam_files)
{
	if (params.barcodes_filename != "" && !std::ifstream(params.barcodes_filename))
		throw std::runtime_error("Can't open barcodes file '" + params.barcodes_filename + "'");

	if (params.config_file_name != "" && !std::ifstream(params.config_file_name))
		throw std::runtime_error("Can't open config file '" + params.config_file_name + "'");

	if (params.genesets_rds != "" && !std::ifstream(params.genesets_rds))
		throw std::runtime_error("Can't open genesets file '" + params.genesets_rds + "'");

	if (params.gtf_filename != "" && !std::ifstream(params.gtf_filename))
		throw std::runtime_error("Can't open GTF file '" + params.gtf_filename + "'");

	if (params.reads_params_names_str != "" && !std::ifstream(params.reads_params_names_str))
		throw std::runtime_error("Can't open reads file '" + params.reads_params_names_str + "'");

	for (auto const &file : bam_files)
	{
		if (!std::ifstream(file))
			throw std::runtime_error("Can't open BAM file '" + file + "'");
	}
}

static void usage()
{
	cerr << "\tindropest: estimate molecular counts per cell" << endl; // Add genes count from config
	cerr << "SYNOPSIS\n";
	cerr <<
	"\tindropest [-t, --text-output] [-m|--merge-cell-tags] [-v|--verbose] [-n | --not-filtered] [-g | --gtf filename]"
	"[-l, --log-prefix logs_name] [-r, --reads-params filename] -c config.xml file1.bam [file2.bam ...] "
	"[-b | --bam-output] [-B | --barcodes filename] [-f, --filled-bam] [-F, --filtered-bam]" << endl;
	cerr << "OPTIONS:\n";
	cerr << "\t-b, --bam-output: print tagged bam files" << endl;
	cerr << "\t-B, --barcodes: path to barcodes file" << endl; //To config
	cerr << "\t-c, --config filename: xml file with estimation parameters" << endl;
	cerr << "\t-f, --filled-bam: bam file already contains genes/barcodes tags" << endl;
	cerr << "\t-F, --filtered-bam: print tagged bam file after the merge and filtration" << endl;
	cerr << "\t-g, --gtf filename: gtf file with genes annotations" << endl;
	cerr << "\t-l, --log-prefix : logs prefix" << endl;
	cerr << "\t-m, --merge-cell-tags : merge linked cell tags" << endl;
	cerr << "\t-n, --not-filtered : print data for all cells" << endl;
	cerr << "\t-o, --output-file filename : output file name" << endl;
	cerr << "\t-p, --parallel number_of_threads : number of threads" << endl;
	cerr << "\t-r, --reads-params filename: file or files with serialized params from tags search step. If there are several files then it should be in quotes and splitted by space" << endl; //DEPRECATED
	cerr << "\t-R, --reads-output: print count matrix for reads and don't use UMI statistics" << endl;
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
			{"barcodes",     	required_argument, 	   0, 'B'},
			{"config",     		required_argument, 0, 'c'},
			{"filled-bam",     	no_argument,       0, 'f'},
			{"filtered-bam",    no_argument,		0, 'F'},
			{"genesets_rds",	required_argument, 0, 'G'},
			{"gtf",     		required_argument, 0, 'g'},
			{"log-prefix",		required_argument, 0, 'l'},
			{"merge-cell-tags", no_argument,       0, 'm'},
			{"not-filtered",	no_argument, 	   0, 'n'},
			{"output-file",     required_argument, 0, 'o'},
			{"parallel",     required_argument, 0, 'p'},
			{"reads-params",     required_argument, 0, 'r'},
			{"reads-output",     no_argument, 		0, 'R'},
			{"text-output",     no_argument,       0, 't'},
			{"verbose",         no_argument,       0, 'v'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "bB:c:fFG:g:l:mno:p:r:Rtv", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'b':
				params.bam_output = true;
				break;
			case 'B':
				params.barcodes_filename = string(optarg);
				break;
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			case 'f' :
				params.filled_bam = true;
				break;
			case 'F' :
				params.filtered_bam_output = true;
				break;
			case 'G' :
				params.genesets_rds = string(optarg);
				break;
			case 'g' :
				params.gtf_filename = string(optarg);
				break;
			case 'l' :
				params.log_prefix = string(optarg);
				break;
			case 'm' :
				params.merge_tags = true;
				break;
			case 'n' :
				params.not_filtered = true;
				break;
			case 'o' :
				params.output_name = string(optarg);
				break;
			case 'p' :
				params.num_of_threads = atoi(optarg);
				break;
			case 'r' :
				params.reads_params_names_str = string(optarg);
				break;
			case 'R' :
				params.reads_output = true;
				break;
			case 't' :
				params.text_output = true;
				break;
			case 'v' :
				params.verbose = true;
				break;
			default:
				cerr << "indropest: unknown arguments passed: '" << c <<"'"  << endl;
				params.cant_parse = true;
				return params;
		}
	}

	if (optind > argc - 1)
	{
		cerr << "indropset: at least one bam file must be supplied" << endl;
		params.cant_parse = true;
	}

	if (params.config_file_name == "")
	{
		cerr << "indropset: config file must be supplied" << endl;
		params.cant_parse = true;
	}

	if (params.filled_bam && params.reads_params_names_str != "")
	{
		cerr << "indropset: only one genes source must be provided (you can't use -r and -f at the same time)" << endl;
		params.cant_parse = true;
	}

	if (!params.merge_tags && params.barcodes_filename != "")
	{
		cerr << "indropset: you should provide barcodes list only for merge procedure (you can't use -B without -m)" << endl;
		params.cant_parse = true;
	}

	if (params.num_of_threads < 1)
	{
		cerr << "The number of threads should be positive" << endl;
		params.cant_parse = true;
	}

	if (params.output_name == "")
	{
		if (params.text_output)
		{
			params.output_name = "cell.counts.txt";
		}
		else
		{
			params.output_name = "cell.counts";
		}
	}

	return params;
}

int main(int argc, char **argv)
{
	std::string command_line;
	for (int i = 0; i < argc; ++i)
	{
		command_line += argv[i];
		command_line += " ";
	}

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

	L_TRACE << command_line;
	vector<string> files;
	while (optind < argc)
	{
		files.push_back(string(argv[optind++]));
	}

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	try
	{
		check_files_existence(params, files);

		Tools::trace_time("Run");
		Estimator estimator(pt.get_child("config.Estimation"));
		CellsDataContainer container = estimator.get_cells_container(files, params.merge_tags, params.bam_output,
		                                                             params.filled_bam, params.reads_params_names_str,
		                                                             params.gtf_filename, params.barcodes_filename);

		if (params.filtered_bam_output)
		{
			BamProcessing::BamController::write_filtered_bam_files(files, params.filled_bam, params.reads_params_names_str,
																   params.gtf_filename, container);
		}

		{
			Results::IndropResult results = estimator.get_results(container, params.not_filtered, params.reads_output);
		
			Tools::trace_time("Done");

			if (params.text_output)
			{
				Results::ResultPrinter::print_text_table(params.output_name, results.cm);
				params.output_name += ".bin";
			}

			results.save_results(container, params.output_name, params.genesets_rds, params.num_of_threads);
		}

		{
			L_TRACE << "Get bad cells results";
			Results::BadCellsStats bad_cells_results = estimator.get_bad_cells_results(container);
			L_TRACE << "Done";
			bad_cells_results.save_rds(container, params.output_name + ".bc.rds");
		}
	}
	catch (std::runtime_error err)
	{
		L_ERR << err.what();
		return 1;
	}
	Tools::trace_time("All done");
}
