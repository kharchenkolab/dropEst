#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Estimation/Estimator.h"
#include "Tools/log_init.h"

using namespace std;
using namespace __gnu_cxx;

struct Params
{
	bool cant_parse = false;
	bool verbose = false;
	bool text_output = false;
	bool merge_tags = false;
	string out_name = "";
	string config_file_name = "";
};

static void usage()
{
	cerr << "\tindropest: estimate molecular counts per cell" << endl;
	cerr << "SYNOPSIS\n";
	cerr <<
	"\tindropest [-g|--min-genes 1000] [-u|--min-umis 10000] [-m|--merge-cell-tags] [-R|--output-r] [-v|--verbose] -c config.xml file1.bam [file2.bam ...]" <<
	endl;
	cerr << "OPTIONS:\n";
	cerr << "\t-o, --output-file filename : output file name" << endl;
	cerr << "\t-t, --text-output : write out text matrix" << endl;
	cerr << "\t-c, --config config.xml : xml file with estimation parameters" << endl;
	cerr << "\t-m, --merge-cell-tags : merge linked cell tags" << endl;
	cerr << "\t-R, --otuput-r : write out RData file" << endl;
}

static Params parse_cmd_params(int argc, char **argv)
{
	Params params;

	int option_index = 0;
	int c;
	static struct option long_options[] = {
			{"verbose",         no_argument,       0, 'v'},
			{"text-output",     no_argument,       0, 't'},
			{"merge-cell-tags", no_argument,       0, 'm'},
			{"output-file",     required_argument, 0, 'o'},
			{"config",     		required_argument, 0, 'c'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "vtmo:c:", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'v' :
				params.verbose = true;
				break;
			case 'm' :
				params.merge_tags = true;
				break;
			case 't' :
				params.text_output = true;
				break;
			case 'o' :
				params.out_name = string(optarg);
				break;
			case 'c' :
				params.config_file_name = string(optarg);
				break;
			default:
				cerr << "indropest: unknown arguments passed" << endl;
				params.cant_parse = true;
				return params;
		}
	}

	if (optind > argc - 1)
	{
		cerr << "indropset: at least one bam file must be supplied" << endl;
		params.cant_parse = true;
		return params;
	}

	if (params.config_file_name == "")
	{
		cerr << "indropset: config file must be supplied" << endl;
		params.cant_parse = true;
		return params;
	}

	if (params.out_name == "")
	{
		if (params.text_output)
		{
			params.out_name = "cell.counts.txt";
		}
		else
		{
			params.out_name = "cell.counts.bin";
		}
	}

	return params;
}

int main(int argc, char **argv)
{
	Params params = parse_cmd_params(argc, argv);

	if (params.cant_parse)
	{
		usage();
		return 1;
	}

	init_log(params.verbose, true, "debug_est.log");

	vector<string> files;
	while (optind < argc)
	{
		files.push_back(string(argv[optind++]));
	}

	boost::property_tree::ptree pt;
	read_xml(params.config_file_name, pt);

	Estimator estimator(files, pt.get_child("config.Estimation"));
	estimator.run(params.merge_tags, params.text_output, params.out_name);
}
