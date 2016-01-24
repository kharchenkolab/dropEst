#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "TagsSearch/SpacerFinder.h"
#include "TagsSearch/TagsFinder.h"
#include "Tools/log_defs.h"
#include "Tools/log_init.h"

using namespace std;
using namespace __gnu_cxx;

struct Params
{
	bool cant_parse = false;
	bool verbose = false;
	string base_name = "";
	string r1_file_name;
	string r2_file_name;
	string config_file_name;
};

static void usage()
{
	cerr << "\tindroptag -- generate tagged indrop fastq files for alignment";
	cerr << "SYNOPSIS";
	cerr << "\tindroptag [-n|--name baseName] [-v|--verbose] read_1.fastq read_2.fastq config.xml";
	cerr << "OPTIONS:";
	cerr << "\t-v, --verbose verbose mode";
	cerr << "\t-n, --name BASE_NAME specify alternative output base name";
}

Params parse_cmd_params(int argc, char **argv)
{
	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"verbose",    no_argument,       0, 'v'},
			{"name",       required_argument, 0, 'n'},
			{0, 0,                            0, 0}
	};

	Params params;
	while ((c = getopt_long(argc, argv, "vn:", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'v' :
				params.verbose = true;
				break;
			case 'n' :
				params.base_name = string(optarg);
				break;
			default:
				cerr << "indroptag: unknown arguments passed";
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

	init_log(params.verbose, true);

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	SpacerFinder spacer_finder(pt.get_child("config.TagsSearch.SpacerSearch"));
	TagsFinder finder(spacer_finder, pt.get_child("config.TagsSearch.TailTrimming"));

	try
	{
		finder.run(params.r1_file_name, params.r2_file_name, params.base_name);
	}
	catch (std::runtime_error &err)
	{
		L_ERR << err.what();
		return 1;
	}


	return 0;
}
