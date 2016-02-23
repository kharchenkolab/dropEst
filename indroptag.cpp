#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "TagsSearch/SpacerFinder.h"
#include "TagsSearch/TagsFinder.h"
#include "Tools/Logs.h"

using namespace std;
using namespace TagsSearch;

struct Params
{
	bool cant_parse;
	bool verbose;
	bool save_reads_names;
	string base_name;
	string log_prefix;
	string r1_file_name;
	string r2_file_name;
	string config_file_name;

	Params() : cant_parse(false), verbose(false), base_name(""), log_prefix("")
	{}
};

static void usage()
{
	cerr << "\tindroptag -- generate tagged indrop fastq files for alignment\n";
	cerr << "SYNOPSIS\n";
	cerr << "\tindroptag [-n|--name baseName] [-v|--verbose] [-s|--save-reads-names] [-l, --log-prefix logs_name] "
					<< "read_1.fastq read_2.fastq config.xml\n";
	cerr << "OPTIONS:\n";
	cerr << "\t-v, --verbose verbose mode\n";
	cerr << "\t-s, --save-reads-names serialize reads parameters to save names\n";
	cerr << "\t-n, --name BASE_NAME specify alternative output base name\n";
	cerr << "\t-l, --log-prefix : logs prefix" << endl;
}

Params parse_cmd_params(int argc, char **argv)
{
	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"verbose",    no_argument,       0, 'v'},
			{"save-reads-names",    no_argument,       0, 's'},
			{"name",       required_argument, 0, 'n'},
			{"log-prefix", required_argument, 0, 'l'},
			{0, 0,                            0, 0}
	};

	Params params;
	while ((c = getopt_long(argc, argv, "vsn:l:", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'v' :
				params.verbose = true;
				break;
			case 's' :
				params.save_reads_names = true;
				break;
			case 'n' :
				params.base_name = string(optarg);
				break;
			case 'l' :
				params.log_prefix = string(optarg);
				break;
			default:
				cerr << "indroptag: unknown arguments passed" << endl;
				usage();
				params.cant_parse = true;
				return params;
		}
	}

	if (optind != argc - 3)
	{
		cerr << "indroptag: two read files and configs file must be provided" << endl;
		usage();
		params.cant_parse = true;
		return params;
	}

	params.r1_file_name = string(argv[optind++]);
	params.r2_file_name = string(argv[optind++]);
	params.config_file_name = string(argv[optind++]);

	if (params.base_name == "")
	{
		params.base_name = params.r2_file_name + ".tagged";
	}

	return params;
}

int main(int argc, char **argv)
{
	Params params = parse_cmd_params(argc, argv);

	if (params.cant_parse)
		return 1;

	if (params.log_prefix.length() != 0)
	{
		params.log_prefix += "_";
	}
	Tools::init_log(params.verbose, false, params.log_prefix + "tag_main.log", params.log_prefix + "tag_debug.log");

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	SpacerFinder spacer_finder(pt.get_child("config.TagsSearch.SpacerSearch"));
	TagsFinder finder(spacer_finder, pt.get_child("config.TagsSearch.TailTrimming"));

	try
	{
		finder.run(params.r1_file_name, params.r2_file_name, params.base_name, params.save_reads_names);
	}
	catch (std::runtime_error &err)
	{
		L_ERR << err.what();
		return 1;
	}


	return 0;
}
