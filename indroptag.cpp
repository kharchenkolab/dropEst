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
	L_ERR << "\tindroptag -- generate tagged indrop fastq files for alignment";
	L_ERR << "SYNOPSIS";
	L_ERR << "\tindroptag [-n|--name baseName] [-v|--verbose] read_1.fastq read_2.fastq config.xml";
	L_ERR << "OPTIONS:";
	L_ERR << "\t-v, --verbose verbose mode";
	L_ERR << "\t-n, --name BASE_NAME specify alternative output base name";
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
	while ((c = getopt_long(argc, argv, "vfn:m:s:l:r:", long_options, &option_index)) != -1)
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
				L_ERR << "indroptag: unknown arguments passed";
				usage();
				params.cant_parse = true;
				return params;
		}
	}

	if (optind != argc - 3)
	{
		L_ERR << "indroptag: two read files and configs file must be provided" << endl;
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
	init_log(boost::log::trivial::debug);

	Params params = parse_cmd_params(argc, argv);

	if (params.cant_parse)
		return 1;

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	SpacerFinder spacer_finder(pt.get_child("config.SpacerSearch"));
	TagsFinder finder(params.verbose, spacer_finder, pt.get_child("config.TailTrimming"));

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
