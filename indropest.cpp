#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>
#include <ctime>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <Estimation/BamProcessing/BamController.h>
#include <Estimation/Merge/MergeStrategyFactory.h>

#include "Estimation/ResultsPrinter.h"
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
	bool merge_tags_precise = false;
	bool reads_output = false;
	bool text_output = false;
	bool quiet = false;
	string config_file_name = "";
	string genes_filename = "";
	string log_prefix = "";
	string output_name = "";
	string reads_params_names_str = ""; //TODO: deprecated. Should be updated after the implementation of quality parsing.
	std::string gene_match_level = CellsDataContainer::Mark::DEFAULT_CODE;
	int max_cells_number = -1;
	int min_genes_after_merge = -1;
	int num_of_threads = 1;
};

static void check_files_existence(const Params &params, const vector<string> &bam_files)
{
	if (params.config_file_name != "" && !std::ifstream(params.config_file_name))
		throw std::runtime_error("Can't open config file '" + params.config_file_name + "'");

	if (params.genes_filename != "" && !std::ifstream(params.genes_filename))
		throw std::runtime_error("Can't open genes file '" + params.genes_filename + "'");

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
	cerr << "\tindropest: estimate molecular counts per cell" << endl;
	cerr << "SYNOPSIS\n";
	cerr <<
	"\tindropest [options] -c config.xml file_1.bam [..., file_n.bam]" << endl;
	cerr << "OPTIONS:\n";
	cerr << "\t-b, --bam-output: print tagged bam files" << endl;
	cerr << "\t-c, --config filename: xml file with estimation parameters" << endl;
	cerr << "\t-C, --cells num: maximal number of output cells" << endl;
	cerr << "\t-f, --filled-bam: bam file already contains genes/barcodes tags" << endl;
	cerr << "\t-F, --filtered-bam: print tagged bam file after the merge and filtration" << endl;
	cerr << "\t-g, --genes filename: file with genes annotations (.bed or .gtf)" << endl;
	cerr << "\t-G, --genes-min num: minimal number of genes in output cells" << endl;
	cerr << "\t-l, --log-prefix : logs prefix" << endl;
	cerr << "\t-L, --gene-match-level :\n"
			"\t\te: count UMIs with exonic reads only;\n"
			"\t\ti: count UMIs with intronic reads only;\n"
			"\t\tE: count UMIs, which have both exonic and not annotated reads;\n"
			"\t\tI: count UMIs, which have both intronic and not annotated reads;\n"
			"\t\tB: count UMIs, which have both exonic and intronic reads;\n"
			"\t\tA: count UMIs, which have exonic, intronic and not annotated reads.\n"
			"\t\tDefault: -L " << Params().gene_match_level << "." << endl;
	cerr << "\t-m, --merge-barcodes : merge linked cell tags" << endl;
	cerr << "\t-M, --merge-barcodes-precise : use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available" << endl;
	cerr << "\t-o, --output-file filename : output file name" << endl;
	cerr << "\t-p, --parallel number_of_threads : number of threads" << endl;
//	cerr << "\t-r, --reads-params filename: file or files with serialized params from tags search step. If there are several files then it should be in quotes and splitted by space" << endl;
	cerr << "\t-R, --reads-output: print count matrix for reads and don't use UMI statistics" << endl;
	cerr << "\t-t, --text-output : write out text matrix" << endl;
	cerr << "\t-q, --quiet : disable logs" << endl;
}

static Params parse_cmd_params(int argc, char **argv)
{
	Params params;

	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"bam-output",     	no_argument, 	   0, 'b'},
			{"config",     		required_argument, 0, 'c'},
			{"cells",     		required_argument, 0, 'C'},
			{"filled-bam",     	no_argument,       0, 'f'},
			{"filtered-bam",    no_argument,		0, 'F'},
			{"genes",     		required_argument, 0, 'g'},
			{"genes-min",     		required_argument, 0, 'G'},
			{"log-prefix",		required_argument, 0, 'l'},
			{"gene-match-level",	required_argument, 0, 'L'},
			{"merge-barcodes",  no_argument,       0, 'm'},
			{"merge-barcodes-precise",  no_argument,       0, 'M'},
			{"not-filtered",	no_argument, 	   0, 'n'},
			{"output-file",     required_argument, 0, 'o'},
			{"parallel",     required_argument, 0, 'p'},
			{"reads-params",     required_argument, 0, 'r'},
			{"reads-output",     no_argument, 		0, 'R'},
			{"text-output",     no_argument,       0, 't'},
			{"quiet",         no_argument,       0, 'q'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "bc:C:fFg:G:l:L:mM:no:p:r:Rtq", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'b':
				params.bam_output = true;
				break;
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			case 'C' :
				params.max_cells_number = atoi(optarg);
				break;
			case 'f' :
				params.filled_bam = true;
				break;
			case 'F' :
				params.filtered_bam_output = true;
				break;
			case 'g' :
				params.genes_filename = string(optarg);
				break;
			case 'G' :
				params.min_genes_after_merge = atoi(optarg);
				break;
			case 'l' :
				params.log_prefix = string(optarg);
				break;
			case 'L' :
				params.gene_match_level = string(optarg);
				break;
			case 'm' :
				params.merge_tags = true;
				break;
			case 'M' :
				params.merge_tags = true;
				params.merge_tags_precise = true;
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
			case 'q' :
				params.quiet = true;
				break;
			default:
				cerr << "indropest: unknown arguments passed: '" << (char)c <<"'"  << endl;
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

	if (params.num_of_threads < 1)
	{
		cerr << "The number of threads should be positive" << endl;
		params.cant_parse = true;
	}

	if (params.genes_filename == "" && params.gene_match_level.find_first_of("eE") == std::string::npos)
	{
		cerr << "indropset: you should provide genes file (-g option) to use intron annotations" << endl;
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
			params.output_name = "cell.counts.rds";
		}
	}

	return params;
}

CellsDataContainer get_cells_container(const vector<string> &files, const Params &params)
{
	boost::property_tree::ptree pt;
	if (!params.config_file_name.empty())
	{
		read_xml(params.config_file_name, pt);
	}

	auto match_levels = CellsDataContainer::Mark::get_by_code(params.gene_match_level);

	Merge::MergeStrategyFactory merge_factory(pt, params.min_genes_after_merge);
	CellsDataContainer container(merge_factory.get_cb_strat(params.merge_tags, params.merge_tags_precise),
	                             merge_factory.get_umi(), match_levels, params.max_cells_number);

	BamProcessing::BamController::parse_bam_files(files, params.bam_output, params.filled_bam,
	                                              params.reads_params_names_str, params.genes_filename, container);

	container.set_initialized();
	container.merge_and_filter();
	return container;
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
	Tools::init_log(!params.quiet, false, params.log_prefix + "est_main.log", params.log_prefix + "est_debug.log");

	L_TRACE << command_line;
	vector<string> files;
	while (optind < argc)
	{
		files.push_back(string(argv[optind++]));
	}

	try
	{
		check_files_existence(params, files);
		Tools::trace_time("Run");
		CellsDataContainer container = get_cells_container(files, params);

		if (params.filtered_bam_output)
		{
			BamProcessing::BamController::write_filtered_bam_files(files, params.filled_bam, params.reads_params_names_str,
																   params.genes_filename, container);
		}

		ResultsPrinter printer(params.text_output, params.reads_output);
		Tools::trace_time("Done");

		if (params.text_output)
		{
			throw std::runtime_error("Text output is not implemented");
		}

		printer.save_results(container, params.output_name);
	}
	catch (std::runtime_error err)
	{
		L_ERR << err.what();
		return 1;
	}
	Tools::trace_time("All done");
}
