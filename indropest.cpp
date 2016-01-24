#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>

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

	int min_genes = 0;
	int min_umis = 0;
	int low_genes = 10; // hard threshold for computational optimizationx
};

static void usage()
{
	cerr << "\tindropest: estimate molecular counts per cell" << endl;
	cerr << "SYNOPSIS\n";
	cerr <<
	"\tindropest [-g|--min-genes 1000] [-u|--min-umis 10000] [-m|--merge-cell-tags] [-R|--output-r] [-v|--verbose] file1.bam [file2.bam ...]" <<
	endl;
	cerr << "OPTIONS:\n";
	cerr << "\t-o, --output-file filename : output file name" << endl;
	cerr << "\t-t, --text-output : write out text matrix" << endl;
	cerr << "\t-g, --min-genes n : output cells with at least n genes" << endl;
	cerr << "\t-u, --min-umis k : output cells with at least k UMIs" << endl;
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
			{"min-genes",       required_argument, 0, 'g'},
			{"min-umis",        required_argument, 0, 'u'},
			{0, 0,                                 0, 0}
	};
	while ((c = getopt_long(argc, argv, "vtmo:g:u:", long_options, &option_index)) != -1)
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
			case 'g' :
				params.min_genes = atoi(optarg);
				break;
			case 'u' :
				params.min_umis = atoi(optarg);
				break;
			case 'o' :
				params.out_name = string(optarg);
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

	if (params.min_genes == 0 && params.min_umis == 0)
	{
		params.min_genes = 1000;
	}

	if (params.min_genes > 0 && params.min_genes < params.low_genes)
	{
		params.low_genes = params.min_genes;
	}
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

	Estimator estimator(files, params.min_genes, params.min_umis, params.low_genes);
	estimator.run(params.merge_tags, params.text_output, params.out_name);
}
