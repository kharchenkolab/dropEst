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
#include "TagsSearch/IndropV1TagsFinder.h"
#include "TagsSearch/TagsFinderBase.h"
#include "TagsSearch/IndropV3TagsFinder.h"
#include <TagsSearch/TextWriter.h>
#include "Tools/Logs.h"
#include <Tools/OmpConcurrentQueue.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

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
	int num_of_threads = 1;
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
	cerr << "\t-p, --parallel number: number of threads\n";
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
	auto const &processing_config = pt.get_child(PROCESSING_CONFIG_PATH, boost::property_tree::ptree());

	if (params.read_files.size() == 4)
	{
		if (params.lib_tag == "")
			throw std::runtime_error("For IndropV3 with library tag, tag (-t option) should be specified");

		return shared_ptr<TagsFinderBase>(
				new IndropV3LibsTagsFinder(params.read_files[0], params.read_files[1], params.read_files[2],
				                           params.read_files[3], params.lib_tag, pt.get_child(BARCODES_CONFIG_PATH),
				                           processing_config, params.save_stats));
	}

	if (params.read_files.size() == 3)
		return shared_ptr<TagsFinderBase>(
				new IndropV3TagsFinder(params.read_files[0], params.read_files[1], params.read_files[2],
				                       pt.get_child(BARCODES_CONFIG_PATH), processing_config, params.save_stats));

	if (params.read_files.size() != 2)
		throw std::runtime_error("Unexpected number of read files: " + std::to_string(params.read_files.size()));

	if (pt.get<std::string>("config.TagsSearch.SpacerSearch.barcode_mask", "") != "")
		return shared_ptr<TagsFinderBase>(
				new FixPosSpacerTagsFinder(params.read_files[0], params.read_files[1], pt.get_child(SPACER_CONFIG_PATH),
				                           processing_config, params.save_stats));

	return shared_ptr<TagsFinderBase>(
			new IndropV1TagsFinder(params.read_files[0], params.read_files[1], pt.get_child(SPACER_CONFIG_PATH),
			                       processing_config, params.save_stats));
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
			{"save-reads-names",    no_argument,       0, 's'},
			{"save-stats",    required_argument,       0, 'S'},
			{"lib-tag",    required_argument, 0, 't'},
			{"quiet",    no_argument,       0, 'q'},
			{0, 0,                            0, 0}
	};

	Params params;
	while ((c = getopt_long(argc, argv, "c:hl:n:p:St:q", long_options, &option_index)) != -1)
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
				params.num_of_threads = atoi(optarg);
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

bool read_bunch(Tools::OmpConcurrentQueue<std::string> &records, shared_ptr<TagsFinderBase> finder,
                size_t number_of_iterations = 10000, size_t records_bunch_size = 5000)
{
	bool ended = false;
	for (size_t i = 0; i < number_of_iterations; ++i)
	{
		if (ended)
			break;

		std::string records_bunch;
		for (size_t record_id = 0; record_id < records_bunch_size; ++record_id)
		{
			FastQReader::FastQRecord record;
			if (finder->file_ended())
			{
				ended = true;
				break;
			}

			if (!finder->get_next_record(record))
			{
				record_id--;
				continue;
			}

			records_bunch += record.to_string();
		}

		if (!records_bunch.empty())
		{
			records.push(records_bunch);
		}
	}

	return ended;
}

void write_all(Tools::OmpConcurrentQueue<std::string> &gzipped, TextWriter &writer)
{
	std::string to_write, cur_text;
	while (gzipped.pop(cur_text))
	{
		to_write += cur_text;
	}
	if (!to_write.empty())
	{
		writer.write(to_write);
	}
}

void gzip_all(Tools::OmpConcurrentQueue<std::string> &input, Tools::OmpConcurrentQueue<std::string> &output, bool bound_out_size)
{
	while (!bound_out_size || !output.full()) // .full shouldn't block the queue. It's just a suggested size.
	{
		std::string records_bunch;

		if (!input.pop(records_bunch))
			break;

		output.push(TextWriter::gzip(records_bunch));
	}
}

void run_parallel(const Params &params, const boost::property_tree::ptree &pt)
{
	shared_ptr<TagsFinderBase> finder = get_tags_finder(params, pt);

	Tools::trace_time("Run");

	bool work_ended = false;
	Tools::OmpConcurrentQueue<std::string> records(100000), gzipped(100000);
	size_t max_file_size = pt
			.get_child(PROCESSING_CONFIG_PATH, boost::property_tree::ptree())
			.get<size_t>("reads_per_out_file", std::numeric_limits<size_t>::max());

	TextWriter writer(params.base_name, "fastq.gz", max_file_size);

	#pragma omp parallel shared(finder, records, gzipped, writer, work_ended)
	while (true)
	{
		bool cur_ended;
		#pragma omp atomic read
		cur_ended = work_ended;

		if (cur_ended && records.empty() && gzipped.empty())
			break;

		// Part 1. Single thread.
		#pragma omp master
		if (!cur_ended)
		{
			cur_ended = read_bunch(records, finder);
			if (cur_ended)
			{
				L_TRACE << finder->results_to_string();
				Tools::trace_time("Reading compleated");
				L_TRACE << "Writing the rest lines";
				#pragma omp atomic write
				work_ended = true;
			}
		}

		#pragma omp atomic read
		cur_ended = work_ended;

		// Part 2. Multithread.
		gzip_all(records, gzipped, !cur_ended);

		#pragma omp master
		write_all(gzipped, writer);
	}

	Tools::trace_time("All done");
	if (params.save_stats)
	{
		save_stats(params.base_name + ".rds", finder);
	}
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

#ifdef _OPENMP
	omp_set_num_threads(params.num_of_threads);
#endif

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
		run_parallel(params, pt);
	}
	catch (std::runtime_error &err)
	{
		time_t ctt = time(0);
		L_ERR << err.what()  << "\nTime: " << asctime(localtime(&ctt));
		return 1;
	}


	return 0;
}
