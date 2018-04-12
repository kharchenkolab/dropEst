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
#include "Tools/UtilFunctions.h"

using namespace std;
using namespace Estimation;
using namespace boost::property_tree;

static const std::string SCRIPT_NAME = "dropest";

struct Params
{
	bool bam_output = false;
	bool cant_parse = false;
	bool filled_bam = false;
	bool filtered_bam_output = false;
	bool merge_tags = false;
	bool merge_tags_precise = false;
	bool pseudoaligner = false;
	bool quiet = false;
	bool reads_output = false;
	bool stats_for_validation = false;
	bool umi_merge = false;
	bool velocyto_matrices = false;
	bool write_matrix = false;
	string config_file_name = "";
	string genes_filename = "";
	string log_prefix = "";
	string output_name = "";
	string read_params_filenames = "";
	std::string gene_match_level = UMI::Mark::DEFAULT_CODE;
	int max_cells_number = -1;
	int min_genes_after_merge = -1;
};

static void check_files_existence(const Params &params, const vector<string> &bam_files)
{
	if (!params.config_file_name.empty() && !std::ifstream(params.config_file_name))
		throw std::runtime_error("Can't open config file '" + params.config_file_name + "'");

	if (!params.genes_filename.empty() && !std::ifstream(params.genes_filename))
		throw std::runtime_error("Can't open genes file '" + params.genes_filename + "'");

	for (auto const &file : bam_files)
	{
		if (!std::ifstream(file))
			throw std::runtime_error("Can't open BAM file '" + file + "'");
	}
}

static void usage()
{
	cerr << SCRIPT_NAME <<": estimate molecular counts per cell\n";
	cerr << "Version: " << VERSION << "\n\n";
	cerr << "SYNOPSIS\n";
	cerr << "\t" << SCRIPT_NAME << " [options] -c config.xml file_1.bam [..., file_n.bam]\n";
	cerr << "OPTIONS:\n";
	cerr << "\t-b, --bam-output: print tagged bam files\n";
	cerr << "\t-c, --config filename: xml file with estimation parameters\n";
	cerr << "\t-C, --cells num: maximal number of output cells\n";
	cerr << "\t-f, --filled-bam: bam file already contains genes/barcodes tags\n";
	cerr << "\t-F, --filtered-bam: print tagged bam file after the merge and filtration\n";
	cerr << "\t-g, --genes filename: file with genes annotations (.bed or .gtf)\n";
	cerr << "\t-G, --genes-min num: minimal number of genes in output cells\n";
	cerr << "\t-h, --help: show this info\n";
	cerr << "\t-l, --log-prefix : logs prefix\n";
	cerr << "\t-L, --gene-match-level :\n"
			"\t\te: count UMIs with exonic reads only;\n"
			"\t\ti: count UMIs with intronic reads only;\n"
			"\t\tE: count UMIs, which have both exonic and not annotated reads;\n"
			"\t\tI: count UMIs, which have both intronic and not annotated reads;\n"
			"\t\tB: count UMIs, which have both exonic and intronic reads;\n"
			"\t\tA: count UMIs, which have exonic, intronic and not annotated reads.\n"
			"\t\tDefault: -L " << Params().gene_match_level << "." << endl;
	cerr << "\t-m, --merge-barcodes : merge linked cell tags" << endl;
	cerr << "\t-M, --merge-barcodes-precise : use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available\n";
	cerr << "\t-o, --output-file filename : output file name\n";
	cerr << "\t-P, --pseudoaligner: use chromosome name as a source of gene id\n";
	cerr << "\t-q, --quiet : disable logs\n";
	cerr << "\t-r, --read-params filenames: file or files with serialized params from tags search step. If there are several files"
	     << ", they should be provided in quotes, separated by space: \"file1.params.gz file2.params.gz file3.params.gz\"" << endl;
	cerr << "\t-R, --reads-output: print count matrix for reads and don't use UMI statistics\n";
	cerr << "\t-u, --merge-umi: apply 'directional' correction of UMI errors. If you want to apply more advanced UMI correction, don’t use ‘-u’, but use follow up R analysis.\n";
	cerr << "\t-V, --velocyto : save separate count matrices for exons, introns and exon/intron spanning reads\n";
	cerr << "\t-w, --write-mtx : write out matrix in MatrixMarket format\n";
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
			{"help",     		no_argument, 0, 'h'},
			{"log-prefix",		required_argument, 0, 'l'},
			{"gene-match-level",	required_argument, 0, 'L'},
			{"merge-barcodes",  no_argument,       0, 'm'},
			{"merge-barcodes-precise",  no_argument,       0, 'M'},
			{"not-filtered",	no_argument, 	   0, 'n'},
			{"output-file",     required_argument, 0, 'o'},
			{"read-params",     required_argument, 0, 'r'},
			{"pseudoaligner",   no_argument, 0, 'P'},
			{"quiet",         no_argument,       0, 'q'},
			{"reads-output",     no_argument, 		0, 'R'},
			{"validation-stats", no_argument,       0, 'S'},
			{"velocyto",     no_argument,       0, 'V'},
			{"write-mtx",     no_argument,       0, 'w'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "bc:C:fFg:G:hl:L:mMno:r:PqRSuVw", long_options, &option_index)) != -1)
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
				params.max_cells_number = int(strtol(optarg, nullptr, 10));
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
				params.min_genes_after_merge = int(strtol(optarg, nullptr, 10));
				break;
			case 'h' :
				usage();
				exit(0);
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
			case 'r' :
				params.read_params_filenames = string(optarg);
				break;
			case 'P' :
				params.pseudoaligner = true;
				break;
			case 'R' :
				params.reads_output = true;
				break;
			case 'q' :
				params.quiet = true;
				break;
			case 'S' :
				params.stats_for_validation = true;
				break;
			case 'V' :
				params.velocyto_matrices = true;
				break;
			case 'u' :
				params.umi_merge = true;
				break;
			case 'w' :
				params.write_matrix = true;
				break;
			default:
				cerr << SCRIPT_NAME << ": unknown arguments passed: '" << (char)c <<"'"  << endl;
				params.cant_parse = true;
				return params;
		}
	}

	if (optind > argc - 1)
	{
		cerr << SCRIPT_NAME << ": at least one bam file must be supplied" << endl;
		params.cant_parse = true;
	}

	if (params.config_file_name.empty())
	{
		cerr << SCRIPT_NAME << ": config file must be supplied" << endl;
		params.cant_parse = true;
	}

	if (params.filled_bam && params.read_params_filenames != "")
	{
		cerr << SCRIPT_NAME << ": only one genes source must be provided (you can't use both -r and -f at the same time)" << endl;
		params.cant_parse = true;
	}

	if (params.genes_filename.empty() && params.gene_match_level.find_first_of("eE") == std::string::npos)
	{
		cerr << SCRIPT_NAME << ": you should provide genes file (-g option) to use intron annotations" << endl;
		params.cant_parse = true;
	}

	if (params.output_name.empty())
	{
		params.output_name = "cell.counts.rds";
	}

	return params;
}

CellsDataContainer get_cells_container(const vector<string> &files, const Params &params,
                                       const ptree &est_config, const BamProcessing::BamController &bam_controller)
{
	auto match_levels = UMI::Mark::get_by_code(params.gene_match_level);

	Merge::MergeStrategyFactory merge_factory(est_config, params.config_file_name, params.min_genes_after_merge);
	CellsDataContainer container(merge_factory.get_cb_strat(params.merge_tags, params.merge_tags_precise),
	                             merge_factory.get_umi(params.umi_merge), match_levels, params.max_cells_number);

	bam_controller.parse_bam_files(files, params.bam_output, container);

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
	Tools::copy_file(params.config_file_name, params.log_prefix + "est_config.dump.xml");

	L_TRACE << command_line;
	L_TRACE << "Version: " << VERSION << ".";
	vector<string> files;

	Tools::init_r();
	while (optind < argc)
	{
		files.emplace_back(argv[optind++]);
	}

	try
	{
		check_files_existence(params, files);
		Tools::trace_time("Run", true);
		ptree estimation_config;
		if (!params.config_file_name.empty())
		{
			ptree pt;
			read_xml(params.config_file_name, pt);
			estimation_config = pt.get_child("config.Estimation", ptree());
		}

		BamProcessing::BamController bam_controller(BamProcessing::BamTags(estimation_config), params.filled_bam,
		                                            params.read_params_filenames, params.genes_filename,
		                                            params.pseudoaligner, estimation_config.get<int>("Other.min_barcode_quality", 0));
		CellsDataContainer container = get_cells_container(files, params, estimation_config, bam_controller);

		if (params.filtered_bam_output)
		{
			bam_controller.write_filtered_bam_files(files, container);
		}

		ResultsPrinter printer(params.write_matrix, params.reads_output, params.stats_for_validation, !params.umi_merge);
		Tools::trace_time("Done");

		printer.save_results(container, params.output_name);

		if (params.velocyto_matrices) {
			printer.save_intron_exon_matrices(container, params.output_name);
		}
	}
	catch (std::runtime_error &err)
	{
		L_ERR << err.what();
		return 1;
	}
	catch (std::logic_error &err)
	{
		L_ERR << err.what();
		return 1;
	}
	catch (std::exception &err)
	{
		L_ERR << err.what();
		return 1;
	}
	Tools::trace_time("All done");
}
