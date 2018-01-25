#include <unordered_set>
#include <string>
#include <vector>

#include <api/BamReader.h>
#include <Tools/Logs.h>
#include <api/BamWriter.h>

using namespace std;
using namespace BamTools;

unordered_set<string> get_read_names(const string &bam_name)
{
	Tools::trace_time("Getting read names from bam: " + bam_name);
	BamReader reader;
	if (!reader.Open(bam_name))
		throw std::runtime_error("Could not open BAM file: " + bam_name);

	BamAlignment alignment;
	unordered_set<string> read_names;

	size_t iter = 0;
	while (reader.GetNextAlignment(alignment))
	{
		if (!alignment.IsMapped() || !alignment.IsPrimaryAlignment())
			continue;

		read_names.emplace(alignment.Name);

		if (++iter % 10000000 == 0)
		{
			L_TRACE << "Iteration " << iter << ": " << read_names.size() << " read names";
		}
	}

	reader.Close();

	L_TRACE << "Done: " << read_names.size() << " read names";
	return read_names;
}

unordered_set<string> get_mixed_reads(const string &org1_bam_name, const string &org2_bam_name)
{
	auto org1_read_names = get_read_names(org1_bam_name);

	Tools::trace_time("Getting read names from bam: " + org2_bam_name);
	BamReader reader;
	if (!reader.Open(org2_bam_name))
		throw std::runtime_error("Could not open BAM file: " + org2_bam_name);

	BamAlignment alignment;
	unordered_set<string> mixed_read_names;

	size_t iter = 0, mixed_reads_num = 0;
	while (reader.GetNextAlignment(alignment))
	{
		if (!alignment.IsMapped() || !alignment.IsPrimaryAlignment())
			continue;

		if (org1_read_names.find(alignment.Name) != org1_read_names.end())
		{
			mixed_read_names.insert(alignment.Name);
			mixed_reads_num++;
		}

		if (++iter % 10000000 == 0)
		{
			L_TRACE << "Iteration " << iter << ": " << mixed_reads_num << " mixed reads";
		}
	}

	reader.Close();

	org1_read_names.clear();

	Tools::trace_time("Done");

	return mixed_read_names;
}

unordered_set<string> get_mixed_reads(const string &bam_name)
{
	L_TRACE << "Reading bam to create list of mixed reads";
	BamReader reader;
	if (!reader.Open(bam_name))
		throw std::runtime_error("Can't open BAM file: " + bam_name);

	BamAlignment alignment;
	unordered_set<string> mouse_reads, human_reads, mixed_reads;

	size_t iter = 0, mixed_reads_num = 0;
	while (reader.GetNextAlignment(alignment))
	{
		if (!alignment.IsMapped() || !alignment.IsPrimaryAlignment())
			continue;

		std::string chr_name = reader.GetReferenceData().at(alignment.RefID).RefName;

		if (++iter % 10000000 == 0)
		{
			L_TRACE << "Iteration " << iter << ": " << mixed_reads_num << " mixed, " << mouse_reads.size()
			        << " mouse, " << human_reads.size() << " human";
		}

		if (mixed_reads.find(alignment.Name) != mixed_reads.end())
		{
			mixed_reads_num++;
			continue;
		}

		if (chr_name[0] == 'h' || chr_name[0] == 'H')
		{
			if (mouse_reads.find(alignment.Name) != mouse_reads.end())
			{
				mixed_reads.emplace(alignment.Name);
				mixed_reads_num++;
			}
			else
			{
				human_reads.emplace(alignment.Name);
			}
		}
		else
		{
			if (human_reads.find(alignment.Name) != human_reads.end())
			{
				mixed_reads.emplace(alignment.Name);
				mixed_reads_num++;
			}
			else
			{
				mouse_reads.insert(alignment.Name);
			}
		}
	}

	reader.Close();

	mouse_reads.clear();
	human_reads.clear();
	L_TRACE << "Done";
	return mixed_reads;
}

void write_filtered_bam(const string &source_bam_name, const string &target_bam_name,
                        const unordered_set<string> &mixed_reads)
{
	Tools::trace_time("Writing filtered bam file", true);

	BamReader reader;
	if (!reader.Open(source_bam_name))
		throw std::runtime_error("Can't open BAM file: " + source_bam_name);

	BamWriter writer;
	if (!writer.Open(target_bam_name, reader.GetConstSamHeader(), reader.GetReferenceData()))
		throw std::runtime_error("Can't open BAM file: " + target_bam_name);

	BamAlignment alignment;
	unordered_set<int32_t> unexpected_chromosome_ids;

	size_t iter = 0, filtered_reads = 0;
	while (reader.GetNextAlignment(alignment))
	{
		if (!alignment.IsMapped() || !alignment.IsPrimaryAlignment())
			continue;

		++iter;
		if (mixed_reads.find(alignment.Name) != mixed_reads.end())
		{
			filtered_reads++;
			continue;
		}

		writer.SaveAlignment(alignment);

		if (iter % 10000000 == 0)
		{
			L_TRACE << "Iteration " << iter << ": " << filtered_reads << " filtered.";
		}
	}

	reader.Close();

	L_TRACE << "Done";
}

int main(int argc, char **argv)
{
	std::string command_line;
	for (int i = 0; i < argc; ++i)
	{
		command_line += argv[i];
		command_line += " ";
	}

	Tools::init_log(true, false, "filt_main.log", "filt_debug.log");
	L_TRACE << command_line;

	vector<string> arguments;

	while (optind < argc)
	{
		arguments.emplace_back(argv[optind++]);
	}

	if (arguments.size() != 1 && arguments.size() != 3)
	{
		L_ERR << "Wrong parameters";
		return 1;
	}

	Tools::trace_time("Run", true);

	try
	{
		unordered_set<string> mixed_reads = arguments.size() == 1
		                                    ? get_mixed_reads(arguments[0])
		                                    : get_mixed_reads(arguments[1], arguments[2]);

		write_filtered_bam(arguments[0], arguments[0] + ".filtered.bam", mixed_reads);
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
	return 0;
}
