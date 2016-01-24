#include "Estimator.h"

#include "IndropResults.h"
#include "Tools/edit_distance.h"
#include "Tools/log_defs.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/range/adaptor/reversed.hpp>

using namespace std;

struct BamRaiiContainer
{
	bam1_t *align_info;
	bamFile bam_input;
	bam_header_t *header;

	BamRaiiContainer(const string &bam_file_name)
	{
		this->align_info = bam_init1();

		if (this->align_info == nullptr)
			throw std::runtime_error("Can't allocate align info");

		this->bam_input = bam_open(bam_file_name.c_str(), "r");
		if (this->bam_input == nullptr)
		{
			bam_destroy1(align_info);
			throw std::runtime_error("Can't open " + bam_file_name);
		}

		this->header = bam_header_read(this->bam_input);
		if (this->header == nullptr)
		{
			bam_close(bam_input);
			bam_destroy1(align_info);
			throw std::runtime_error("Can't read header");
		}
	}

	~BamRaiiContainer()
	{
		bam_header_destroy(header);
		bam_close(bam_input);
		bam_destroy1(align_info);
	}
};

Estimator::Estimator(const std::vector<std::string> &files, int min_genes, int min_umis, int low_genes)
		: min_genes(min_genes)
		, min_umis(min_umis)
		, low_genes(low_genes)
{
	for (auto const &name : files)
	{
		this->parse_bam_file(name);
	}
}

void Estimator::parse_bam_file(const string &bam_file_name)
{
	L_TRACE << "reading " << bam_file_name;

	BamRaiiContainer c(bam_file_name);

	long total_reads = 0, exonic_reads = 0;
	while (bam_read1(c.bam_input, c.align_info) >= 0)
	{
		total_reads++;
		if (total_reads % 1000000 == 0)
		{
			L_TRACE << "Total " << total_reads << " reads processed";
		}

		string chr(c.header->target_name[c.align_info->core.tid]);
		uint8_t *ptr = bam_aux_get(c.align_info, "GE");
		if (ptr == nullptr)
		{ // classify non-exonic read
			nonexone_chrs[chr]++;
			continue;
		}

		exonic_reads++;
		// parse out the tags
		string qname(bam1_qname(c.align_info));
		size_t umi_start_pos = qname.rfind('#');
		if (umi_start_pos == string::npos)
		{
			L_ERR << "WARNING: unable to parse out UMI in qname: " << qname;
			continue;
		}
		string umi = qname.substr(umi_start_pos + 1);
		size_t cell_tag_start_pos = qname.rfind('!', umi_start_pos);
		if (cell_tag_start_pos == string::npos)
		{
			L_ERR << "WARNING: unable to parse out cell tag in qname: " << qname;
			continue;
		}
		string cell_tag = qname.substr(cell_tag_start_pos + 1, umi_start_pos - cell_tag_start_pos - 1);

		if (c.align_info->core.l_qseq < this->read_prefix_length)
		{
			L_ERR << "WARNING: read is shorter than read_prefix_length. total_reads=" << total_reads << endl;
			continue;
		}

		string gene(bam_aux2Z(ptr));

		L_DEBUG << qname << " cell:" << cell_tag << " UMI:" << umi << " prefix:"
				<< Estimator::get_iseq_verbose(c.align_info, this->read_prefix_length) << "\tXF:" << gene;

		// insert into the cb map
		auto res = cb_ids.emplace(cell_tag, this->cb_genes.size());
		if (res.second)
		{ // new cb
			this->cb_genes.push_back(SHHM());
			this->cb_names.push_back(cell_tag);
		}
		int cb_id = res.first->second;
		this->cb_genes[cb_id][gene][umi]++;
		string umig = umi + gene; // +iseq
		this->umig_cbs[umig][cb_id]++;
		//string umip=umi+prefix; // +iseq
		//umip_cbs[umip][cb_id]++;

		exone_chrs[chr]++;

		L_DEBUG << "CB/UMI=" << this->cb_genes[cb_id][gene][umi] << " gene=" << this->cb_genes[cb_id][gene].size()
				<< " CB=" << this->cb_genes[cb_id].size() << " UMIg=" << this->umig_cbs[umig].size();

//		if (total_reads > 1000000)
//			break;
	}

	L_TRACE << " done (" << total_reads << " total reads; " << exonic_reads << " exonic reads; "
			<< this->cb_genes.size() << " cell barcodes)";
}

void Estimator::run(bool merge_tags, bool text_output, const string &output_name)
{
	// count the number of UMIs per cb
	ints_t merge_n;
	ints_t unmerged_cbs;
	i_counter_t cb_genes_count_by_id(this->count_cb_genes()); // <ngenes,cb_id> pairs

	if (merge_tags)
	{
		this->merge_genes_with_tags(cb_genes_count_by_id, merge_n, unmerged_cbs);
		cb_genes_count_by_id = this->count_cb_genes(false);
	}
	else
	{
		this->merge_genes_without_tags(cb_genes_count_by_id, unmerged_cbs);
	}

	this->umig_cbs.clear(); // free up some memory

	L_TRACE << unmerged_cbs.size() << " valid (with >=" << min_genes << " genes) cells with ";

	s_counter_t gene_counts(this->count_genes(unmerged_cbs));

	L_TRACE << "compiling count matrix ... ";
	names_t cell_names(this->get_unmerged_names(unmerged_cbs));

	names_t gene_names;
	ints_t umis;
	this->parse_genes(unmerged_cbs, gene_counts, gene_names, umis);

	L_TRACE << "Done";

	string bin_output_name = output_name;
	if (text_output)
	{
		bin_output_name += ".bin";
		this->print_text_output(output_name, gene_names, cell_names, umis);
	}

	CountMatrix cm(cell_names, gene_names, umis);

	IndropResult reults = this->get_indrop_results(cm, cb_genes_count_by_id, unmerged_cbs, merge_n);
	this->print_bin_output(bin_output_name, reults);
}

Estimator::i_counter_t Estimator::count_cb_genes(bool logs) const
{
	i_counter_t cb_genes_count_by_id; // <ngenes,cb_id> pairs
	for (int i = 0; i < this->cb_genes.size(); i++)
	{
		size_t genes_count = this->cb_genes[i].size();
		if (genes_count >= this->low_genes)
		{
			cb_genes_count_by_id.push_back(make_pair(genes_count, i));
		}
	}

	if (logs)
	{
		L_TRACE << cb_genes_count_by_id.size() << " CBs with more than " << this->low_genes << " genes";
	}

	sort(cb_genes_count_by_id.begin(), cb_genes_count_by_id.end(),
		 [](const pair<int, int> &p1, const pair<int, int> &p2) { return p1.first < p2.first; });

	if (logs)
	{
		L_TRACE << this->get_cb_count_top_verbose(cb_genes_count_by_id);
	}

	return cb_genes_count_by_id;
}

string Estimator::get_iseq_verbose(bam1_t *align_info, size_t read_prefix_length) const
{
	// unpack first part of the read
	string iseq(read_prefix_length, ' ');

	uint8_t *s = bam1_seq(align_info);

	for (int i = 0; i < read_prefix_length; i++)
	{
		iseq[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
	}

	return iseq;
}

string Estimator::get_cb_top_verbose(const ints_t &unmerged_cbs) const
{
	stringstream ss;
	if (unmerged_cbs.size() > 0)
	{
		ss << "top CBs:\n";
		for (int i = 0; i < min(unmerged_cbs.size(), this->top_print_size); i++)
		{
			ss << this->cb_genes[unmerged_cbs[i]].size() << "\t" << this->cb_names[unmerged_cbs[i]] << "\n";
		}
	}
	else
	{
		cout << "no valid CBs found\n";
	}

	return ss.str();
}

string Estimator::get_cb_count_top_verbose(const i_counter_t &cb_genes_count_by_id) const
{
	stringstream ss;
	if (cb_genes_count_by_id.size() > 0)
	{
		ss << "top CBs:\n";
		size_t low_border = cb_genes_count_by_id.size() - min(cb_genes_count_by_id.size(), this->top_print_size);
		for (size_t i = cb_genes_count_by_id.size() - 1; i > low_border; --i)
		{
			ss << cb_genes_count_by_id[i].first << "\t" << this->cb_names[cb_genes_count_by_id[i].second] << "\n";
		}
	}
	else
	{
		ss << "no valid CBs found\n";
	}

	return ss.str();
}

std::string Estimator::get_genes_top_verbose(const Estimator::s_counter_t &genes) const
{
	ostringstream ss;
	ss << "top genes:\n";
	for (size_t i = 0; i < min(genes.size(), this->top_print_size); i++)
	{
		ss << genes[i].first << '\t' << genes[i].second << "\n";
	}
	return ss.str();
}

Estimator::s_counter_t Estimator::count_genes(const ints_t &unmerged_cbs) const
{
	SIHM counter;
	for (int cb: unmerged_cbs)
	{
		for (auto const &gm: this->cb_genes[cb])
		{
			counter[gm.first] += gm.second.size();
		}
	}

	L_TRACE << counter.size() << " genes";

	s_counter_t gene_counts(counter.begin(), counter.end());
	sort(gene_counts.begin(), gene_counts.end(),
		 [](const pair<string, int> &p1, const pair<string, int> &p2) { return p1.second > p2.second; });

	L_TRACE << this->get_genes_top_verbose(gene_counts);
	return gene_counts;
}

void Estimator::merge_genes_with_tags(const i_counter_t &cb_genes_count_by_id, ints_t &merge_n, ints_t &unmerged_cbs)
{
	// cb merging
	int merges_count = 0;
	// reassigned barcodes ids
	ints_t cb_reassigned(this->cb_genes.size());
	std::iota(cb_reassigned.begin(), cb_reassigned.end(), 0);

	ISIHM cb_reassigned_to;

	L_TRACE << "merging linked tags ";

	int tag_index = 0;
	for (auto const &gene_count : cb_genes_count_by_id)
	{ // iterate through the minimally-slected CBs, from low to high counts
		tag_index++;
		if (tag_index % 1000 == 0)
		{
			L_TRACE << "Total " << tag_index << " tags merged";
		}
		int gene_id = gene_count.second;

		IIHM umig_top;
		int umig_n = 0;
		for (auto const &i: this->cb_genes[gene_id])
		{
			string gene = i.first;
			//cout<<gene<<":[";
			for (auto j: i.second)
			{
				//cout<<j.first<<":(";
				//look up umig
				string umig = j.first + gene;
				for (auto const &k: this->umig_cbs[umig])
				{
					//cout<<k.first<<" ";
					if (cb_reassigned[k.first] != gene_id && this->cb_genes[cb_reassigned[k.first]].size() > gene_count.first)
					{
						umig_top[cb_reassigned[k.first]]++;
					}
				}
				umig_n++;
				//cout<<") ";
			}
			//cout<<"]"<<endl;
		}
		// get top umig, its edit distance
		int top_cb = -1;
		double top_cb_fraction = -1;
		int top_cb_genes_count = -1;
		for (auto l: umig_top)
		{
			double cb_f = (((double) l.second) / ((double) umig_n));
			if (cb_f > top_cb_fraction || (cb_f == top_cb_fraction &&
										this->cb_genes[l.first].size() > top_cb_genes_count))
			{
				top_cb = l.first;
				top_cb_fraction = cb_f;
				top_cb_genes_count = this->cb_genes[l.first].size();
			}
		}
		bool merged = false;
		if (top_cb > 0)
		{
			//cout<<"\ttop cb: "<<this->cb_names[top_cb]<<" ("<<top_cb_genes_count<<" genes)="<<top_cb_fraction<<" ";
			// check if the top candidate is valid for merging
			if (top_cb_fraction > this->min_merge_fraction)
			{
				int ed = edit_distance(this->cb_names[top_cb].c_str(), this->cb_names[gene_id].c_str());
				//cout<<" ed="<<ed<<endl;
				if (ed < this->max_merge_edit_distance)
				{
					// do the merge
					merged = true;
					//cout<<"merging "<<kid<<" ("<<this->cb_genes[kid].size()<<" genes) into "<<top_cb<<" ("<<top_cb_genes_count<<" genes) ";
					merge_n.push_back(-1 * gene_count.first);
					merges_count++;
					// merge the actual data
					for (auto l: this->cb_genes[gene_id])
					{
						for (auto m: l.second)
						{
							this->cb_genes[top_cb][l.first][m.first] += m.second;
						}
					}
					// establish new mapping
					cb_reassigned[gene_id] = top_cb; // set reassignment mapping
					cb_reassigned_to[top_cb].insert(gene_id); // reassign current cell

					// transfer mapping of the cbs previously mapped to kid
					auto k = cb_reassigned_to.find(gene_id);
					if (k != cb_reassigned_to.end())
					{
						for (auto m: k->second)
						{
							cb_reassigned[m] = top_cb; // update reassignment mapping
						}
						cb_reassigned_to[top_cb].insert(k->first); // reassign to the new cell
					}
					//cout<<this->cb_genes[top_cb].size()<<" genes"<<endl;
					//cout<<"\t"<<top_cb<<": "<<this->cb_genes[top_cb].size()<<" genes, "<<cb_reassigned_to[top_cb].size()<<" reassigned cbs"<<endl;
				}
			}
		}
		if (!merged)
		{
			//cout<<" not merging"<<endl;
			merge_n.push_back(gene_count.first);
			if (gene_count.first >= min_genes)
			{ // only record cells that are passing min_genes threshold
				unmerged_cbs.push_back(gene_id);
			}
		}
	}

	L_INFO << " done (" << merges_count << " merges performed)" << endl;

	if (unmerged_cbs.size() > 1)
	{
		reverse(unmerged_cbs.begin(), unmerged_cbs.end());
	}

	L_TRACE << this->get_cb_top_verbose(unmerged_cbs);
}

void Estimator::merge_genes_without_tags(const Estimator::i_counter_t &cb_genes_count_by_id,
										 Estimator::ints_t &unmerged_cbs) const
{
	// just pick out all sufficiently informative cells
	for (auto const &gene_count : boost::adaptors::reverse(cb_genes_count_by_id))
	{
		if (gene_count.first < min_genes) // only record cells that are passing min_genes threshold
			break;

		unmerged_cbs.push_back(gene_count.second);
	}
}

Estimator::names_t Estimator::get_unmerged_names(const Estimator::ints_t &unmerged_cbs) const
{
	names_t cell_names;
	cell_names.reserve(unmerged_cbs.size());
	for (int i: unmerged_cbs)
	{
		cell_names.push_back(this->cb_names[i]);
	}
	return cell_names;
}

void Estimator::parse_genes(const Estimator::ints_t &unmerged_cbs, const Estimator::s_counter_t &gene_counts,
							Estimator::names_t &gene_names, Estimator::ints_t &umis) const
{
	unsigned long size = unmerged_cbs.size() * gene_counts.size();
	umis.resize(size);
	gene_names.reserve(gene_counts.size());

	for (int i = 0; i < gene_counts.size(); i++)
	{
		const string gene_name = gene_counts[i].first;
		gene_names.push_back(gene_name);
		for (int j = 0; j < unmerged_cbs.size(); j++)
		{
			auto const &cur_genes = this->cb_genes[unmerged_cbs[j]];
			auto res = cur_genes.find(gene_name);

			umis[(i * unmerged_cbs.size()) + j] = res == cur_genes.end() ? 0 : res->second.size();
		}
	}
}

void Estimator::print_text_output(const std::string &output_name, const Estimator::names_t &gene_names,
								  const Estimator::names_t &cell_names, const Estimator::ints_t &umis) const
{
	L_TRACE << "Writing output matrix to " << output_name << " ";

	ofstream output_file(output_name.c_str(), ios_base::out);
	// header
	output_file << "gene";
	for (const string &cell_name: cell_names)
	{
		output_file << '\t' << cell_name;
	}
	output_file << endl;
	for (int i = 0; i < gene_names.size(); i++)
	{
		output_file << gene_names[i];
		for (int j = 0; j < cell_names.size(); j++)
		{
			int umi = umis[(i * cell_names.size()) + j];
			output_file << '\t' << umi;
		}
		output_file << endl;
	}
	//output_file.pop();
	//ofs.close();
	output_file.close();

	L_TRACE << "Done";
}

IndropResult Estimator::get_indrop_results(const CountMatrix cm, const i_counter_t &cb_genes_count_by_id,
										   const ints_t &unmerged_cbs, const ints_t &merge_n) const
{
	L_TRACE << "compiling diagnostic stats: ";

	// calculate average reads/UMI / cell
	vector<double> reads_per_umis(this->get_reads_per_umis(unmerged_cbs));
	L_TRACE << "reads/UMI";

	ints_t umig_coverage(this->get_umig_coverage(cb_genes_count_by_id));
	L_TRACE << "UMIg coverage";

	// serialize non-exonic chromosome counts
	names_t non_exon_count_names;
	ints_t non_exon_counts;
	this->split_pairs(this->nonexone_chrs, non_exon_count_names, non_exon_counts);

	names_t exon_count_names;
	ints_t exon_counts;
	this->split_pairs(this->exone_chrs, exon_count_names, exon_counts);

	return IndropResult(cm, non_exon_counts, non_exon_count_names, reads_per_umis, umig_coverage,
						 exon_counts, exon_count_names, merge_n);
}

std::vector<double> Estimator::get_reads_per_umis(const ints_t &unmerged_cbs) const
{
	vector<double> reads_per_umis;
	for (int j = 0; j < unmerged_cbs.size(); j++)
	{
		int umis_count = 0;
		double reads_pre_umi = 0.0;
		for (auto generec: this->cb_genes[unmerged_cbs[j]])
		{
			for (auto umirec: generec.second)
			{
				reads_pre_umi += umirec.second;
				umis_count++;
			}
		}
		reads_pre_umi /= umis_count;
		reads_per_umis.push_back(reads_pre_umi);
	}

	return reads_per_umis;
}

Estimator::ints_t Estimator::get_umig_coverage(const i_counter_t &cb_genes_count_by_id) const
{
	ints_t umig_coverage;
	unordered_set<string> umigs_seen;
	for (auto const &gene_count : boost::adaptors::reverse(cb_genes_count_by_id))
	{
		int newumigs = 0;
		for (auto generec: this->cb_genes[gene_count.second])
		{
			for (auto umirec: generec.second)
			{
				string umig = umirec.first + generec.first;
				auto res = umigs_seen.emplace(umig);
				if (res.second)
				{
					newumigs++;
				}
			}
		}
		umig_coverage.push_back(newumigs);
	}
	return umig_coverage;
}

void Estimator::split_pairs(const SIHM &base, names_t &out1, ints_t out_2) const
{
	for (auto i:base)
	{
		out1.push_back(i.first);
		out_2.push_back(i.second);
	}
}

void Estimator::print_bin_output(const string &bin_output_name, const IndropResult &results) const
{
	L_TRACE << "Done";
	L_TRACE << "writing binary results to " << bin_output_name;

	ofstream bin_out_file(bin_output_name.c_str(), ios_base::out | ios_base::binary);
	boost::archive::binary_oarchive oa(bin_out_file);
	oa << results;
	bin_out_file.close();

	L_TRACE << "All done";
}
