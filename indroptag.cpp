#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem/path.hpp>
#include <TagsSearch/IndropV3LibsTagsFinder.h>

#include "TagsSearch/FixPosSpacerTagsFinder.h"
#include "TagsSearch/SpacerTagsFinder.h"
#include "TagsSearch/TagsFinderAbstract.h"
#include "TagsSearch/IndropV3TagsFinder.h"

#include "Tools/Logs.h"

using namespace std;
using namespace TagsSearch;

struct Params
{
	bool cant_parse = false;
	bool verbose = false;
	bool save_reads_names = false;
	string base_name = "";
	string log_prefix = "";
	string lib_tag = "";
	string config_file_name = "";
	vector<string> read_files = vector<string>();
};

static void usage()
{
	cerr << "\tindroptag -- generate tagged indrop fastq files for alignment\n";
	cerr << "SYNOPSIS\n";
	cerr << "\tindroptag [-n|--name baseName] [-v|--verbose] [-s|--save-reads-names] [-l, --log-prefix logs_name] "
					<< "-c config.xml read_1.fastq read_2.fastq [read_3.fastq] [library_tags.fastq]\n";
	cerr << "OPTIONS:\n";
	cerr << "\t-c, --config filename: xml file with indroptag parameters" << endl;
	cerr << "\t-l, --log-prefix : logs prefix" << endl;
	cerr << "\t-n, --name BASE_NAME specify alternative output base name\n";
	cerr << "\t-s, --save-reads-names serialize reads parameters to save names\n";
	cerr << "\t-t, --lib-tag library tag (for IndropV3 with library tag only)\n";
	cerr << "\t-v, --verbose verbose mode\n";
}

shared_ptr<TagsFinderAbstract> get_tags_finder(const Params &params, const boost::property_tree::ptree &pt)
{
	auto files_processor = make_shared<FilesProcessor>(params.read_files, params.base_name, params.save_reads_names);

	if (params.read_files.size() == 4)
	{
		if (params.lib_tag == "")
			throw std::runtime_error("For IndropV3 with library tag, tag (-t option) should be specified");

		return shared_ptr<TagsFinderAbstract>(new IndropV3LibsTagsFinder(files_processor, params.lib_tag, 2, //TODO to parameters
		                                                                 pt.get_child("config.TagsSearch.BarcodesSearch"),
		                                                                 pt.get_child("config.TagsSearch.TailTrimming")));
	}

	if (params.read_files.size() == 3)
		return shared_ptr<TagsFinderAbstract>(new IndropV3TagsFinder(files_processor,
																	pt.get_child("config.TagsSearch.BarcodesSearch"),
																	pt.get_child("config.TagsSearch.TailTrimming")));

	if (params.read_files.size() != 2)
		throw std::runtime_error("Unexpected number of read files: " + std::to_string(params.read_files.size()));

	if (pt.get<std::string>("config.TagsSearch.SpacerSearch.barcode_mask", "") != "")
		return shared_ptr<TagsFinderAbstract>(new FixPosSpacerTagsFinder(files_processor,
																	 pt.get_child("config.TagsSearch.SpacerSearch"),
																	 pt.get_child("config.TagsSearch.TailTrimming")));

	return shared_ptr<TagsFinderAbstract>(new SpacerTagsFinder(files_processor,
														   pt.get_child("config.TagsSearch.SpacerSearch"),
														   pt.get_child("config.TagsSearch.TailTrimming")));
}


Params parse_cmd_params(int argc, char **argv)
{
	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"config",     required_argument, 0, 'c'},
			{"log-prefix", required_argument, 0, 'l'},
			{"name",       required_argument, 0, 'n'},
			{"save-reads-names",    no_argument,       0, 's'},
			{"lib-tag",    required_argument, 0, 't'},
			{"verbose",    no_argument,       0, 'v'},
			{0, 0,                            0, 0}
	};

	Params params;
	while ((c = getopt_long(argc, argv, "c:l:n:sv", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			case 'l' :
				params.log_prefix = string(optarg);
				break;
			case 'n' :
				params.base_name = string(optarg);
				break;
			case 's' :
				params.save_reads_names = true;
				break;
			case 't' :
				params.lib_tag= string(optarg);
				break;
			case 'v' :
				params.verbose = true;
				break;
			default:
				cerr << "indroptag: unknown arguments passed" << endl;
				usage();
				params.cant_parse = true;
				return params;
		}
	}

	if (params.config_file_name == "")
	{
		cerr << "indroptag: config file must be supplied" << endl;
		params.cant_parse = true;
	}

	if (optind == argc)
	{
		cerr << "indroptag: read files must be supplied" << endl;
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

	if (params.log_prefix.length() != 0)
	{
		params.log_prefix += "_";
	}
	Tools::init_log(params.verbose, false, params.log_prefix + "tag_main.log", params.log_prefix + "tag_debug.log");

	L_TRACE << command_line;

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	time_t ctt = time(0);
	try
	{
		shared_ptr<TagsFinderAbstract> finder = get_tags_finder(params, pt);

		L_TRACE << "Run: " << asctime(localtime(&ctt));
		finder->run(params.save_reads_names);
		ctt = time(0);
		L_TRACE << "All done: " << asctime(localtime(&ctt));
	}
	catch (std::runtime_error &err)
	{
		ctt = time(0);
		L_ERR << err.what()  << "\nTime: " << asctime(localtime(&ctt));
		return 1;
	}


	return 0;
}