#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem/path.hpp>
#include <TagsSearch/IndropV3LibsTagsFinder.h>
#include <RInside.h>
#include <Tools/UtilFunctions.h>

#include "TagsSearch/FixPosSpacerTagsFinder.h"
#include "TagsSearch/SpacerTagsFinder.h"
#include "TagsSearch/TagsFinderBase.h"
#include "TagsSearch/IndropV3TagsFinder.h"

#include "Tools/Logs.h"

using namespace std;
using namespace TagsSearch;

const std::string PROCESSING_CONFIG_PATH = "config.TagsSearch.Processing";
const std::string BARCODES_CONFIG_PATH = "config.TagsSearch.BarcodesSearch";
const std::string SPACER_CONFIG_PATH = "config.TagsSearch.SpacerSearch";

static const std::string SCRIPT_NAME = "droptag";

struct Params
{
	bool cant_parse = false;
	bool save_reads_names = false;
	bool quiet = false;
	bool save_stats = false;
	string base_name = "";
	string config_file_name = "";
	string log_prefix = "";
	string lib_tag = "";
	vector<string> read_files = vector<string>();
};

void save_stats(const string &out_filename, shared_ptr<TagsFinderBase> tags_finder);

static void usage()
{
	cerr << SCRIPT_NAME << " -- generate tagged fastq files for alignment\n";
	cerr << "Version: " << VERSION << "\n\n";
	cerr << "SYNOPSIS\n";
	cerr << "\t" << SCRIPT_NAME << " [options] "
	     << "-c config.xml barcode_reads.fastq [barcode_umi_reads.fastq] gene_reads.fastq [library_tags.fastq]\n";
	cerr << "OPTIONS:\n";
	cerr << "\t-c, --config filename: xml file with " << SCRIPT_NAME << " parameters\n";
	cerr << "\t-h, --help: show this info\n";
	cerr << "\t-l, --log-prefix prefix: logs prefix\n";
	cerr << "\t-n, --name name: alternative output base name\n";
//	cerr << "\t-s, --save-reads-names : serialize reads parameters to save names\n";
	cerr << "\t-S, --save-stats : save stats to rds file\n";
	cerr << "\t-t, --lib-tag library tag : (for IndropV3 with library tag only)\n";
	cerr << "\t-q, --quiet : disable logs\n";
}

static void check_files_existence(const Params &params)
{
	if (params.config_file_name != "" && !std::ifstream(params.config_file_name))
		throw std::runtime_error("Can't open config file '" + params.config_file_name + "'");

	for (auto const &file : params.read_files)
	{
		if (!std::ifstream(file))
			throw std::runtime_error("Can't open reads file '" + file + "'");
	}
}

shared_ptr<TagsFinderBase> get_tags_finder(const Params &params, const boost::property_tree::ptree &pt)
{
	auto files_processor = make_shared<FilesProcessor>(params.read_files, params.base_name, params.save_reads_names);

	auto const &processing_config = pt.get_child(PROCESSING_CONFIG_PATH, boost::property_tree::ptree());

	if (params.read_files.size() == 4)
	{
		if (params.lib_tag == "")
			throw std::runtime_error("For IndropV3 with library tag, tag (-t option) should be specified");

		return shared_ptr<TagsFinderBase>(new IndropV3LibsTagsFinder(files_processor, params.lib_tag,
	                                                                 pt.get_child(BARCODES_CONFIG_PATH),
		                                                             processing_config));
	}

	if (params.read_files.size() == 3)
		return shared_ptr<TagsFinderBase>(new IndropV3TagsFinder(files_processor, pt.get_child(BARCODES_CONFIG_PATH),
		                                                         processing_config));

	if (params.read_files.size() != 2)
		throw std::runtime_error("Unexpected number of read files: " + std::to_string(params.read_files.size()));

	if (pt.get<std::string>("config.TagsSearch.SpacerSearch.barcode_mask", "") != "")
		return shared_ptr<TagsFinderBase>(new FixPosSpacerTagsFinder(files_processor, pt.get_child(SPACER_CONFIG_PATH),
		                                                             processing_config));

	return shared_ptr<TagsFinderBase>(new SpacerTagsFinder(files_processor, pt.get_child(SPACER_CONFIG_PATH),
	                                                       processing_config));
}


Params parse_cmd_params(int argc, char **argv)
{
	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"config",     required_argument, 0, 'c'},
			{"help", no_argument, 0, 'h'},
			{"log-prefix", required_argument, 0, 'l'},
			{"name",       required_argument, 0, 'n'},
			{"save-reads-names",    no_argument,       0, 's'},
			{"save-stats",    required_argument,       0, 'S'},
			{"lib-tag",    required_argument, 0, 't'},
			{"quiet",    no_argument,       0, 'q'},
			{0, 0,                            0, 0}
	};

	Params params;
	while ((c = getopt_long(argc, argv, "c:hl:n:St:q", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			case 'l' :
				params.log_prefix = string(optarg);
				break;
			case 'h' :
				usage();
				exit(0);
			case 'n' :
				params.base_name = string(optarg);
				break;
			case 's' :
				params.save_reads_names = true;
				break;
			case 'S' :
				params.save_stats = true;
				break;
			case 't' :
				params.lib_tag= string(optarg);
				break;
			case 'q' :
				params.quiet = true;
				break;
			default:
				cerr << SCRIPT_NAME << ": unknown arguments passed" << endl;
				usage();
				params.cant_parse = true;
				return params;
		}
	}

	if (params.config_file_name == "")
	{
		cerr << SCRIPT_NAME << ": config file must be supplied" << endl;
		params.cant_parse = true;
	}

	if (optind == argc)
	{
		cerr << SCRIPT_NAME << ": read files must be supplied" << endl;
		usage();
		params.cant_parse = true;
		return params;
	}

	while (optind != argc)
	{
		params.read_files.push_back(string(argv[optind++]));
	}

	if (params.base_name == "")
	{
		params.base_name = boost::filesystem::path(params.read_files.back()).filename().generic_string() + ".tagged";
	}

	return params;
}

void save_stats(const string &out_filename, shared_ptr<TagsFinderBase> tags_finder)
{
	using namespace Rcpp;

	Tools::trace_time("Writing R data to " + out_filename);
	RInside *R = Tools::init_r();

	(*R)["d"] = List::create(
		_["reads_per_cb"] = wrap(tags_finder->num_reads_per_cb())
	);

	R->parseEvalQ("saveRDS(d, '" + out_filename + "')");
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
		return 1;

	check_files_existence(params);

	if (params.log_prefix.length() != 0)
	{
		params.log_prefix += "_";
	}
	Tools::init_log(!params.quiet, false, params.log_prefix + "tag_main.log", params.log_prefix + "tag_debug.log");

	L_TRACE << command_line;
	L_TRACE << "Version: " << VERSION << ".";

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	try
	{
		shared_ptr<TagsFinderBase> finder = get_tags_finder(params, pt);

		Tools::trace_time("Run");
		finder->run(params.save_reads_names, params.save_stats);
		Tools::trace_time("All done");
		if (params.save_stats)
		{
			save_stats(params.base_name + ".rds", finder);
		}
	}
	catch (std::runtime_error &err)
	{
		time_t ctt = time(0);
		L_ERR << err.what()  << "\nTime: " << asctime(localtime(&ctt));
		return 1;
	}


	return 0;
}