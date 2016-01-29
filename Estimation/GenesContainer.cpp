#include "GenesContainer.h"

#include "Tools/edit_distance.h"
#include "Tools/log_defs.h"

#include <unordered_set>

#include <boost/range/adaptor/reversed.hpp>

using namespace std;

const double EPS = 0.00001;

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

GenesContainer::GenesContainer(const names_t &files, bool merge_tags, size_t read_prefix_length,
							   double min_merge_fraction, int min_genes_before_merge, int min_genes_after_merge,
							   int max_merge_edit_distance, size_t top_print_size)
	: top_print_size(top_print_size)
{
	s_i_hash_t cb_ids;
	s_ii_hash_t umig_cells_counts;
	for (auto const &name : files)
	{
		this->parse_bam_file(name, read_prefix_length, cb_ids, umig_cells_counts);
	}

	this->_cells_genes_counts_sorted = this->count_cells_genes(min_genes_before_merge);

	if (merge_tags)
	{
		this->merge_genes(umig_cells_counts, min_merge_fraction, min_genes_after_merge, max_merge_edit_distance);

		this->_cells_genes_counts_sorted = this->count_cells_genes(min_genes_before_merge, false);
	}
	else
	{
		for (auto const &gene_count : boost::adaptors::reverse(this->_cells_genes_counts_sorted))
		{
			if (gene_count.count < min_genes_after_merge)
				break;

			this->_filtered_cells.push_back(gene_count.index);
		}
	}
}

void GenesContainer::merge_genes(const s_ii_hash_t &umig_cells_counts, double min_merge_fraction,
								 int min_genes_after_merge, int max_merge_edit_distance)
{
	typedef std::unordered_map<int, unordered_set<int> > ISIHM;
	// cb merging
	int merges_count = 0;
	// reassigned barcodes ids
	ints_t cb_reassigned(this->_cells_genes.size());
	iota(cb_reassigned.begin(), cb_reassigned.end(), 0);

	ISIHM cb_reassigned_to;

	L_TRACE << "merging linked tags ";

	int tag_index = 0;
	for (auto const &gene_count : this->_cells_genes_counts_sorted)
	{ // iterate through the minimally-selected CBs, from low to high counts
		tag_index++;
		if (tag_index % 1000 == 0)
		{
			L_TRACE << "Total " << tag_index << " tags merged";
		}

		i_i_hash_t umig_top;
		size_t umigs_count = this->get_umig_top(cb_reassigned, gene_count, umig_cells_counts, umig_top);
		size_t cell_id = gene_count.index;

		// get top umig, its edit distance
		int top_cell_ind = -1;
		double top_cell_fraction = -1;
		int top_cell_genes_count = -1;
		for (auto const &cell: umig_top)
		{
			double cb_fraction = (((double) cell.second) / ((double) umigs_count));
			if (cb_fraction - top_cell_fraction > EPS || (abs(cb_fraction - top_cell_fraction) < EPS &&
														  this->_cells_genes[cell.first].size() > top_cell_genes_count))
			{
				top_cell_ind = cell.first;
				top_cell_fraction = cb_fraction;
				top_cell_genes_count = this->_cells_genes[cell.first].size();
			}
		}

		bool merged = false;
		if (top_cell_ind > 0)
		{
			// check if the top candidate is valid for merging
			if (top_cell_fraction > min_merge_fraction)
			{
				int ed = Tools::edit_distance(this->_genes_names[top_cell_ind].c_str(), this->_genes_names[cell_id].c_str());
				if (ed < max_merge_edit_distance)
				{
					// do the merge
					merged = true;

					this->_stats.add_merge_count(-1 * gene_count.count);
					merges_count++;

					// merge the actual data
					for (auto const &gene: this->_cells_genes[cell_id])
					{
						for (auto const &umi_count: gene.second)
						{
							this->_cells_genes[top_cell_ind][gene.first][umi_count.first] += umi_count.second;
						}
					}
					// establish new mapping
					cb_reassigned[cell_id] = top_cell_ind; // set reassignment mapping
					cb_reassigned_to[top_cell_ind].insert(cell_id); // reassign current cell

					// transfer mapping of the cbs previously mapped to kid
					auto k = cb_reassigned_to.find(cell_id);
					if (k != cb_reassigned_to.end())
					{
						for (auto m: k->second)
						{
							cb_reassigned[m] = top_cell_ind; // update reassignment mapping
						}
						cb_reassigned_to[top_cell_ind].insert(k->first); // reassign to the new cell
					}
				}
			}
		}
		if (!merged)
		{
			this->_stats.add_merge_count(gene_count.count);
			if (gene_count.count >= min_genes_after_merge)
			{
				this->_filtered_cells.push_back(cell_id);
			}
		}
	}

	if (this->_filtered_cells.size() > 1)
	{
		reverse(this->_filtered_cells.begin(), this->_filtered_cells.end());
	}

	L_INFO << "Done (" << merges_count << " merges performed)" << endl;
}

void GenesContainer::parse_bam_file(const string &bam_file_name, size_t read_prefix_length,
									s_i_hash_t &cells_ids, s_ii_hash_t &umig_cells_counts)
{
	L_TRACE << "Reading " << bam_file_name;

	BamRaiiContainer c(bam_file_name);

	long total_reads = 0, exonic_reads = 0;
	while (bam_read1(c.bam_input, c.align_info) >= 0)
	{
		total_reads++;
		if (total_reads % 1000000 == 0)
		{
			L_TRACE << "Total " << total_reads << " reads processed";
		}

		if (c.align_info->core.l_qseq < read_prefix_length)
		{
			L_ERR << "WARNING: read is shorter than read_prefix_length. total_reads=" << total_reads << endl;
			continue;
		}

		string chr_name(c.header->target_name[c.align_info->core.tid]);
		uint8_t *ptr = bam_aux_get(c.align_info, "GE");
		if (ptr == nullptr)
		{
			this->_stats.inc_nonexone_chr_reads(chr_name);
			continue;
		}

		string qname, cell_tag, umi;
		this->parse_read_name(c.align_info, qname, cell_tag, umi);

		string gene(bam_aux2Z(ptr));

		L_DEBUG << qname << " cell:" << cell_tag << " UMI:" << umi << " prefix:"
		<< GenesContainer::get_iseq_verbose(c.align_info, read_prefix_length) << "\tXF:" << gene;

		auto res = cells_ids.emplace(cell_tag, this->_cells_genes.size());
		if (res.second)
		{ // new cb
			this->_cells_genes.push_back(genes_t());
			this->_genes_names.push_back(cell_tag);
		}
		int cell_id = res.first->second;
		this->_cells_genes[cell_id][gene][umi]++;
		string umig = umi + gene; // +iseq
		umig_cells_counts[umig][cell_id]++;

		exonic_reads++;
		this->_stats.inc_exone_chr_reads(chr_name);

		L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi] << " gene=" << this->_cells_genes[cell_id][gene].size()
				<< " CB=" << this->_cells_genes[cell_id].size() << " UMIg=" << umig_cells_counts[umig].size();
	}

	L_TRACE << "Done (" << total_reads << " total reads; " << exonic_reads << " exonic reads; "
			<< this->_cells_genes.size() << " cell barcodes)";
}

GenesContainer::i_counter_t GenesContainer::count_cells_genes(int min_genes_before_merge, bool logs) const
{
	i_counter_t cells_genes_counts; // <genes_count,cell_id> pairs
	for (size_t i = 0; i < this->_cells_genes.size(); i++)
	{
		size_t genes_count = this->_cells_genes[i].size();
		if (genes_count >= min_genes_before_merge)
		{
			cells_genes_counts.push_back(IndexedCount{i, genes_count});
		}
	}

	if (logs)
	{
		L_TRACE << cells_genes_counts.size() << " CBs with more than " << min_genes_before_merge << " genes";
	}

	sort(cells_genes_counts.begin(), cells_genes_counts.end(),
		 [](const IndexedCount &p1, const IndexedCount &p2) { return p1.count < p2.count; });

	if (logs)
	{
		L_TRACE << this->get_cb_count_top_verbose(cells_genes_counts);
	}

	return cells_genes_counts;
}

bool GenesContainer::parse_read_name(const bam1_t *align_info, string &read_name, string &cell_tag, string &umi) const
{
	read_name = bam1_qname(align_info);
	size_t umi_start_pos = read_name.rfind('#');
	if (umi_start_pos == string::npos)
	{
		L_ERR << "WARNING: unable to parse out UMI in read_name: " << read_name;
		return false;
	}

	size_t cell_tag_start_pos = read_name.rfind('!', umi_start_pos);
	if (cell_tag_start_pos == string::npos)
	{
		L_ERR << "WARNING: unable to parse out cell tag in read_name: " << read_name;
		return false;
	}

	umi = read_name.substr(umi_start_pos + 1);
	cell_tag = read_name.substr(cell_tag_start_pos + 1, umi_start_pos - cell_tag_start_pos - 1);

	return true;
}


string GenesContainer::get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const
{
	stringstream ss;
	if (cells_genes_counts.size() > 0)
	{
		ss << "top CBs:\n";
		size_t low_border = cells_genes_counts.size() - min(cells_genes_counts.size(), this->top_print_size);
		for (size_t i = cells_genes_counts.size() - 1; i > low_border; --i)
		{
			ss << cells_genes_counts[i].count << "\t" << this->_genes_names[cells_genes_counts[i].index] << "\n";
		}
	}
	else
	{
		ss << "no valid CBs found\n";
	}

	return ss.str();
}

string GenesContainer::get_iseq_verbose(bam1_t *align_info, size_t read_prefix_length) const
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


const Stats &GenesContainer::stats() const
{
	return this->_stats;
}

const GenesContainer::genes_t &GenesContainer::cell_genes(size_t index) const
{
	return this->_cells_genes[index];
}

const GenesContainer::i_counter_t &GenesContainer::cells_genes_counts_sorted() const
{
	return this->_cells_genes_counts_sorted;
}

const std::string &GenesContainer::gene_name(size_t index) const
{
	return this->_genes_names[index];
}

const GenesContainer::ids_t &GenesContainer::filtered_cells() const
{
	return this->_filtered_cells;
}

size_t GenesContainer::get_umig_top(const ints_t &cb_reassigned, const IndexedCount &gene_count,
									const s_ii_hash_t &umig_cells_counts, i_i_hash_t &umig_top) const
{
	size_t umigs_count = 0;
	for (auto const &gene: this->_cells_genes[gene_count.index])
	{
		string gene_name = gene.first;

		for (auto const &umi_count: gene.second)
		{
			string umig = umi_count.first + gene_name;
			for (auto const &cells_counts: umig_cells_counts.at(umig))
			{
				if (cb_reassigned[cells_counts.first] == gene_count.index)
					continue;

				if (this->_cells_genes[cb_reassigned[cells_counts.first]].size() > gene_count.count)
				{
					umig_top[cb_reassigned[cells_counts.first]]++;
				}
			}
			umigs_count++;
		}
	}
	return umigs_count;
}
