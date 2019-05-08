#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem/path.hpp>
#include <TagsSearch/IndropV3LibsTagsFinder.h>
#include <Tools/UtilFunctions.h>

#include "TagsSearch/FixPosSpacerTagsFinder.h"
#include "TagsSearch/IndropV1TagsFinder.h"
#include "TagsSearch/TagsFinderBase.h"
#include "TagsSearch/IndropV3TagsFinder.h"
#include <TagsSearch/ConcurrentGzWriter.h>
#include <TagsSearch/IClipTagsFinder.h>
#include <TagsSearch/MultipleBarcodesTagsFinder.h>
#include <Rcpp.h>
#include "Tools/Logs.h"

using namespace std;
using namespace TagsSearch;
using namespace boost::property_tree;

const std::string CONFIG_PATH = "config.TagsSearch";
const std::string PROCESSING_CONFIG_PATH = CONFIG_PATH + ".Processing";
const std::string MULTIPLE_BARCODES_CONFIG_PATH = CONFIG_PATH + ".MultipleBarcodeSearch";
const std::string BARCODES_CONFIG_PATH = CONFIG_PATH + ".BarcodesSearch";
const std::string SPACER_CONFIG_PATH = CONFIG_PATH + ".SpacerSearch";

static const std::string SCRIPT_NAME = "droptag";

struct Params
{
	bool cant_parse = false;
	bool save_reads_params = false;
	bool quiet = false;
	bool save_stats = false;
	int num_of_threads = 1;
        int reads_per_out_file = -1;
	string base_name = "";
	string config_file_name = "";
	string log_prefix = "";
	string lib_tag = "";
	vector<string> read_files = vector<string>();
};

void save_stats(const string &out_filename, const shared_ptr<TagsFinderBase> &tags_finder);

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
	cerr << "\t-p, --parallel number: number of threads\n";
	cerr << "\t-s, --save-reads-params : serialize reads parameters to save quality info\n";
	cerr << "\t-S, --save-stats : save stats to rds file\n";
	cerr << "\t-r, --reads-per-out-file : maximum number of reads per output file; (0 - unlimited). Overrides corresponding xml parameter.\n";
	cerr << "\t-t, --lib-tag library tag : (for IndropV3 with library tag only)\n";
	cerr << "\t-q, --quiet : disable logs\n";
}

static void check_files_existence(const Params &params)
{
	if (!params.config_file_name.empty() && !std::ifstream(params.config_file_name))
		throw std::runtime_error("Can't open config file '" + params.config_file_name + "'");

	for (auto const &file : params.read_files)
	{
		if (!std::ifstream(file))
			throw std::runtime_error("Can't open file with reads: '" + file + "'");
	}
}

/// Returns instance of TagsFinderBase based on the config file and number of input files, initializing it according to CLI parameters
/// @param params CLI parameters
/// @param pt specific node in the config file (config/TagsSearch)
shared_ptr<TagsFinderBase> get_tags_finder(const Params &params, const ptree &pt)
{
	auto const &config = pt.get_child(CONFIG_PATH);
	std::string protocol_type = config.get<std::string>("protocol", "");
	if (protocol_type.empty())
		throw std::runtime_error("Protocol is empty. Please, specify it in the config (TagsSearch/protocol)");

	auto const &processing_config = pt.get_child(PROCESSING_CONFIG_PATH, ptree());

	size_t max_records_per_file = params.reads_per_out_file;
	if(max_records_per_file==-1) max_records_per_file = processing_config.get<size_t>("reads_per_out_file", 0);
	auto writer = std::make_shared<ConcurrentGzWriter>(params.base_name, "fastq.gz", max_records_per_file);

	const std::string input_files_num_error_text = "Unexpected number of read files: " +
			std::to_string(params.read_files.size()) +
			" for protocol \"" + protocol_type + "\"";

	if (protocol_type == "indrop3")
	{
		if (params.read_files.size() == 4)
		{
			if (params.lib_tag.empty())
				throw std::runtime_error("For IndropV3 with library tag, tag (-t option) should be specified");

			return shared_ptr<TagsFinderBase>(
					new IndropV3LibsTagsFinder(params.read_files, params.lib_tag, pt.get_child(BARCODES_CONFIG_PATH),
					                           processing_config, writer, params.save_stats, params.save_reads_params));
		}

		if (params.read_files.size() != 3)
			throw std::runtime_error(input_files_num_error_text);

		return shared_ptr<TagsFinderBase>(
				new IndropV3TagsFinder(params.read_files, pt.get_child(BARCODES_CONFIG_PATH), processing_config,
				                       writer, params.save_stats, params.save_reads_params));
	}

	if (protocol_type == "10x") // 10x has the same format of files as indrop v3
	{
		if (params.read_files.size() != 3)
			throw std::runtime_error(input_files_num_error_text);

		return shared_ptr<TagsFinderBase>(
				new IndropV3TagsFinder(params.read_files, pt.get_child(BARCODES_CONFIG_PATH), processing_config,
				                       writer, params.save_stats, params.save_reads_params));
	}

	if (protocol_type == "indrop")
	{
		if (params.read_files.size() != 2)
			throw std::runtime_error(input_files_num_error_text);

		if (!pt.get<std::string>(SPACER_CONFIG_PATH + ".barcode_mask", "").empty())
			return shared_ptr<TagsFinderBase>(
					new FixPosSpacerTagsFinder(params.read_files, pt.get_child(SPACER_CONFIG_PATH), processing_config,
					                           writer, params.save_stats, params.save_reads_params));

		return shared_ptr<TagsFinderBase>(
				new IndropV1TagsFinder(params.read_files, pt.get_child(SPACER_CONFIG_PATH), processing_config,
				                       writer, params.save_stats, params.save_reads_params));
	}

	if (protocol_type == "iclip")
	{
		if (params.read_files.size() != 1)
			throw std::runtime_error(input_files_num_error_text);

		return shared_ptr<TagsFinderBase>(
				new IClipTagsFinder(params.read_files, pt.get_child(BARCODES_CONFIG_PATH), processing_config,
				                    writer, params.save_stats, params.save_reads_params));
	}

	if (protocol_type == "split_seq")
	{
		if (params.read_files.size() != 2)
			throw std::runtime_error(input_files_num_error_text);

		return shared_ptr<TagsFinderBase>(
				new MultipleBarcodesTagsFinder(params.read_files, pt.get_child(MULTIPLE_BARCODES_CONFIG_PATH), processing_config,
						writer, params.save_stats, params.save_reads_params));
	}

    throw std::runtime_error("Unknown protocol: '" + protocol_type + "'. Please, check it in the config (TagsSearch/protocol)");
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
			{"parallel",   required_argument, 0, 'p'},
			{"reads-per-out-file",   required_argument, 0, 'r'},
			{"save-reads-params",    no_argument,       0, 's'},
			{"save-stats",    required_argument,       0, 'S'},
			{"lib-tag",    required_argument, 0, 't'},
			{"quiet",    no_argument,       0, 'q'},
			{0, 0,                            0, 0}
	};

	Params params;
	while ((c = getopt_long(argc, argv, "c:hl:n:p:r:sSt:q", long_options, &option_index)) != -1)
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
			case 'p' :
				params.num_of_threads = int(strtol(optarg, nullptr, 10));
				break;
			case 's' :
				params.save_reads_params = true;
				break;
			case 'r' :
				params.reads_per_out_file = int(strtol(optarg, nullptr, 10));  
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

	if (params.config_file_name.empty())
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
		params.read_files.emplace_back(argv[optind++]);
	}

	if (params.base_name.empty())
	{
		params.base_name = boost::filesystem::path(params.read_files.back()).filename().generic_string() + ".tagged";
	}

	return params;
}

void save_stats(const string &out_filename, const shared_ptr<TagsFinderBase> &tags_finder)
{
	using namespace Rcpp;

	Tools::trace_time("Writing R data to " + out_filename);
	auto R = Tools::init_r();

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
	Tools::copy_file(params.config_file_name, params.log_prefix + "tag_config.dump.xml");

	L_TRACE << command_line;
	L_TRACE << "Version: " << VERSION << ".";

	ptree pt;
	xml_parser::read_xml(params.config_file_name, pt);

	try
	{
		shared_ptr<TagsFinderBase> finder = get_tags_finder(params, pt);

		Tools::trace_time("Run");

		// All work is done here
		finder->run(params.num_of_threads);

		Tools::trace_time("All done");
		if (params.save_stats)
		{
			save_stats(params.base_name + ".rds", finder);
		}
	}
	catch (std::runtime_error &err)
	{
		time_t ctt = time(nullptr);
		L_ERR << err.what()  << "\nTime: " << asctime(localtime(&ctt));
		return 1;
	}


	return 0;
}
