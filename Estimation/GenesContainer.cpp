#include "GenesContainer.h"

#include "Tools/UtilFunctions.h"
#include "Tools/Logs.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <fstream>

//#include <boost/range/adaptor/reversed.hpp>

using namespace std;

namespace Estimation
{
	const double EPS = 0.00001;

	struct BamRaiiContainer
	{
		bam1_t *align_info;
		bamFile bam_input;
		bam_header_t *header;

		BamRaiiContainer(const string &bam_file_name)
		{
			this->align_info = bam_init1();

			if (this->align_info == 0)
				throw runtime_error("Can't allocate align info");

			this->bam_input = bam_open(bam_file_name.c_str(), "r");
			if (this->bam_input == 0)
			{
				bam_destroy1(align_info);
				throw runtime_error("Can't open " + bam_file_name);
			}

			this->header = bam_header_read(this->bam_input);
			if (this->header == 0)
			{
				bam_close(bam_input);
				bam_destroy1(align_info);
				throw runtime_error("Can't read header");
			}
		}

		~BamRaiiContainer()
		{
			bam_header_destroy(header);
			bam_close(bam_input);
			bam_destroy1(align_info);
		}
	};

	GenesContainer::GenesContainer(size_t read_prefix_length, double min_merge_fraction, int min_genes_before_merge,
								   int min_genes_after_merge, int max_merge_edit_distance, size_t top_print_size,
								   const std::string &reads_params_name)
		: _top_print_size(top_print_size)
		, _min_merge_fraction(min_merge_fraction)
		, _min_genes_after_merge(min_genes_after_merge)
		, _min_genes_before_merge(min_genes_before_merge)
		, _max_merge_edit_distance(max_merge_edit_distance)
		, _read_prefix_length(read_prefix_length)
		, _use_names_map(reads_params_name.length() != 0)
	{
		if (this->_use_names_map)
		{
			std::ifstream ifs(reads_params_name);
			boost::archive::binary_iarchive ia(ifs);
			ia >> this->_reads_params;
		}
	}

	void GenesContainer::init(const std::vector<std::string> &files, bool merge_tags)
	{
		s_i_hash_t cb_ids;
		s_ii_hash_t umig_cells_counts;
//	for (auto const &name : files)
		for (names_t::const_iterator name_it = files.begin(); name_it != files.end(); ++name_it)
		{
			this->parse_bam_file(*name_it, cb_ids, umig_cells_counts);
		}

		this->_cells_genes_counts_sorted = this->count_cells_genes();

		if (merge_tags)
		{
			this->merge_genes(umig_cells_counts);

			this->_cells_genes_counts_sorted = this->count_cells_genes(false);
		}
		else
		{
//		for (auto const &gene_count : boost::adaptors::reverse(this->_cells_genes_counts_sorted))
			for (i_counter_t::const_reverse_iterator gene_count_it = this->_cells_genes_counts_sorted.rbegin();
				 gene_count_it != this->_cells_genes_counts_sorted.rend(); ++gene_count_it)
			{
				if (gene_count_it->count < this->_min_genes_after_merge)
					break;

				this->_filtered_cells.push_back(gene_count_it->index);
			}
		}
	}

	void GenesContainer::merge_genes(const s_ii_hash_t &umig_cells_counts)
	{
		// cb merging
		int merges_count = 0;
		// reassigned barcodes ids
		ints_t cb_reassigned(this->_cells_genes.size());
		for (int i = 0; i < cb_reassigned.size(); ++i)
		{
			cb_reassigned[i] = i;
		}

		ISIHM cb_reassigned_to;

		L_TRACE << "merging linked tags ";

		int tag_index = 0;
//	for (auto const &gene_count : this->_cells_genes_counts_sorted)
		for (i_counter_t::const_iterator gene_count_it = this->_cells_genes_counts_sorted.begin();
			 gene_count_it != this->_cells_genes_counts_sorted.end(); ++gene_count_it)
		{ // iterate through the minimally-selected CBs, from low to high counts
			const IndexedCount &cur_gene = *gene_count_it;
			tag_index++;
			if (tag_index % 1000 == 0)
			{
				L_TRACE << "Total " << tag_index << " tags merged";
			}

			i_i_hash_t umig_top;
			size_t umigs_count = this->get_umig_top(cb_reassigned, cur_gene, umig_cells_counts, umig_top);

			// get top umig, its edit distance
			int top_cell_ind = -1;
			double top_cell_fraction = -1;
			long top_cell_genes_count = -1;
//		for (auto const &cell: umig_top)
			for (i_i_hash_t::const_iterator cell_it = umig_top.begin(); cell_it != umig_top.end(); ++cell_it)
			{
				const pair<int, int> &cell = *cell_it;
				double cb_fraction = cell.second / (double) umigs_count;
				if (cb_fraction - top_cell_fraction > EPS || (abs(cb_fraction - top_cell_fraction) < EPS &&
															  this->_cells_genes[cell.first].size() >
															  top_cell_genes_count))
				{
					top_cell_ind = cell.first;
					top_cell_fraction = cb_fraction;
					top_cell_genes_count = this->_cells_genes[cell.first].size();
				}
			}

			if (this->merge(top_cell_ind, top_cell_fraction, cur_gene, cb_reassigned, cb_reassigned_to))
			{
				merges_count++;
			}
			else
			{
				this->_stats.add_merge_count(cur_gene.count);
				if (cur_gene.count >= this->_min_genes_after_merge)
				{
					this->_filtered_cells.push_back(cur_gene.index);
				}
			}
		}

		this->_stats.merge(cb_reassigned, this->_cells_names);

		if (this->_filtered_cells.size() > 1)
		{
			reverse(this->_filtered_cells.begin(), this->_filtered_cells.end());
		}

		L_INFO << "Done (" << merges_count << " merges performed)" << endl;
	}

	bool GenesContainer::merge(int top_cell_ind, double top_cell_fraction, const IndexedCount &gene_count,
							   ints_t &cb_reassigned, ISIHM &cb_reassigned_to)
	{
		if (top_cell_ind < 0)
			return false;

		// check if the top candidate is valid for merging
		if (top_cell_fraction < this->_min_merge_fraction)
			return false;

		size_t cell_id = gene_count.index;
		int ed = Tools::edit_distance(this->_cells_names[top_cell_ind].c_str(), this->_cells_names[cell_id].c_str());
		if (ed >= this->_max_merge_edit_distance)
			return false;

		// do the merge
		this->_stats.add_merge_count(-1 * gene_count.count);

		// merge the actual data
//	for (auto const &gene: this->_cells_genes[cell_id])
		genes_t &cell_genes = this->_cells_genes[cell_id];
		for (genes_t::const_iterator gene_it = cell_genes.begin(); gene_it != cell_genes.end(); ++gene_it)
		{
//		for (auto const &umi_count: gene_it->second)
			for (s_i_hash_t::const_iterator umi_count_it = gene_it->second.begin();
				 umi_count_it != gene_it->second.end(); ++umi_count_it)
			{
				this->_cells_genes[top_cell_ind][gene_it->first][umi_count_it->first] += umi_count_it->second;
			}
		}

		this->reassign(cell_id, top_cell_ind, cb_reassigned, cb_reassigned_to);

		return true;
	}

	void GenesContainer::reassign(size_t cell_id, int target_cell_id, ints_t &cb_reassigned,
								  ISIHM &cb_reassigned_to) const
	{
		cb_reassigned[cell_id] = target_cell_id; // set reassignment mapping
		cb_reassigned_to[target_cell_id].insert(cell_id); // reassign current cell

		// transfer mapping of the cbs previously mapped to kid
		ISIHM::iterator k = cb_reassigned_to.find(cell_id);
		if (k == cb_reassigned_to.end())
			return;

//	for (auto m: k->second)
		for (i_set_t::const_iterator it = k->second.begin(); it != k->second.end(); ++it)
		{
			cb_reassigned[*it] = target_cell_id; // update reassignment mapping
			cb_reassigned_to[target_cell_id].insert(*it);
		}

		k->second.clear();
		cb_reassigned_to[target_cell_id].insert(k->first); // reassign to the new cell
	}

	void GenesContainer::parse_bam_file(const string &bam_file_name, s_i_hash_t &cells_ids,
										s_ii_hash_t &umig_cells_counts)
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

			if (c.align_info->core.l_qseq < this->_read_prefix_length)
			{
				L_ERR << "WARNING: read is shorter than read_prefix_length. total_reads=" << total_reads << endl;
				continue;
			}

			string chr_name(c.header->target_name[c.align_info->core.tid]);

			string qname, cell_barcode, umi;
			if (!this->parse_read_name(c.align_info, qname, cell_barcode, umi))
				continue;

			uint8_t *ptr = bam_aux_get(c.align_info, "GE");
			if (ptr == 0)
			{
				this->_stats.inc_cell_chr_umi(chr_name, cell_barcode, Stats::NON_EXONE);
				continue;
			}

			string gene(bam_aux2Z(ptr));

			L_DEBUG << qname << " cell:" << cell_barcode << " UMI:" << umi << " prefix:"
					<< GenesContainer::get_iseq_verbose(c.align_info) << "\tXF:" << gene;

			pair<s_i_hash_t::iterator, bool> res = cells_ids.emplace(cell_barcode, this->_cells_genes.size());
			if (res.second)
			{ // new cb
				this->_cells_genes.push_back(genes_t());
				this->_cells_names.push_back(cell_barcode);
			}
			int cell_id = res.first->second;
			this->_cells_genes[cell_id][gene][umi]++;
			string umig = umi + gene; // +iseq
			umig_cells_counts[umig][cell_id]++;

			exonic_reads++;
			this->_stats.inc(Stats::READS_BY_UMIG, cell_barcode + "_" + umig);
			this->_stats.inc_cell_chr_umi(chr_name, cell_barcode, Stats::EXONE);

			L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi] << " gene=" <<
					this->_cells_genes[cell_id][gene].size()
					<< " CB=" << this->_cells_genes[cell_id].size() << " UMIg=" << umig_cells_counts[umig].size();
		}

		L_TRACE << "Done (" << total_reads << " total reads; " << exonic_reads << " exonic reads; "
				<< this->_cells_genes.size() << " cell barcodes)";
	}

	GenesContainer::i_counter_t GenesContainer::count_cells_genes(bool logs) const
	{
		i_counter_t cells_genes_counts; // <genes_count,cell_id> pairs
		for (size_t i = 0; i < this->_cells_genes.size(); i++)
		{
			size_t genes_count = this->_cells_genes[i].size();
			if (genes_count >= this->_min_genes_before_merge)
			{
				cells_genes_counts.push_back(IndexedCount(i, genes_count));
			}
		}

		if (logs)
		{
			L_TRACE << cells_genes_counts.size() << " CBs with more than " << this->_min_genes_before_merge << " genes";
		}

		sort(cells_genes_counts.begin(), cells_genes_counts.end(), IndexedCount::counts_comp);

		if (logs)
		{
			L_TRACE << this->get_cb_count_top_verbose(cells_genes_counts);
		}

		return cells_genes_counts;
	}

	bool GenesContainer::parse_read_name(const bam1_t *align_info, string &read_name, string &cell_barcode,
										 string &umi) const
	{
		read_name = bam1_qname(align_info);
		Tools::ReadsParameters params;

		if (this->_use_names_map)
		{
			auto iter = this->_reads_params.find(read_name);
			if (iter == this->_reads_params.end())
			{
				L_ERR << "WARNING: can't find read_name in map: " << read_name;
				return false;
			}

			params = iter->second;
			if (params.is_empty()) {
				L_ERR << "WARNING: empty parameters for read_name: " << read_name;
				return false;
			}
		}
		else
		{
			try
			{
				params = Tools::ReadsParameters(read_name);
			}
			catch (std::runtime_error &error)
			{
				L_ERR << error.what();
				return false;
			}
		}

		cell_barcode = params.cell_barcode();
		umi = params.umi_barcode();
		return true;
	}


	string GenesContainer::get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const
	{
		stringstream ss;
		if (cells_genes_counts.size() > 0)
		{
			ss << "top CBs:\n";
			size_t low_border = cells_genes_counts.size() - min(cells_genes_counts.size(), this->_top_print_size);
			for (size_t i = cells_genes_counts.size() - 1; i > low_border; --i)
			{
				ss << cells_genes_counts[i].count << "\t" << this->_cells_names[cells_genes_counts[i].index] << "\n";
			}
		}
		else
		{
			ss << "no valid CBs found\n";
		}

		return ss.str();
	}

	string GenesContainer::get_iseq_verbose(bam1_t *align_info) const
	{
		// unpack first part of the read
		string iseq(this->_read_prefix_length, ' ');

		uint8_t *s = bam1_seq(align_info);

		for (int i = 0; i < this->_read_prefix_length; i++)
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

	const string &GenesContainer::cell_name(size_t index) const
	{
		return this->_cells_names[index];
	}

	const GenesContainer::ids_t &GenesContainer::filtered_cells() const
	{
		return this->_filtered_cells;
	}

	size_t GenesContainer::get_umig_top(const ints_t &cb_reassigned, const IndexedCount &cur_gene,
										const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const
	{
		size_t umigs_count = 0;
//	for (auto const &gene: this->_cells_genes[cur_gene.index])
		const genes_t &cur_cell_genes = this->_cells_genes[cur_gene.index];
		for (genes_t::const_iterator gene_it = cur_cell_genes.begin(); gene_it != cur_cell_genes.end(); ++gene_it)
		{
			const string &gene_name = gene_it->first;
			const s_i_hash_t &umis = gene_it->second;

//		for (auto const &umi_count: umis)
			for (s_i_hash_t::const_iterator umi_it = umis.begin(); umi_it != umis.end(); ++umi_it)
			{
				string umig = umi_it->first + gene_name;
				const i_i_hash_t &umig_cells = umigs_cells_counts.at(umig);
				for (i_i_hash_t::const_iterator cells_it = umig_cells.begin(); cells_it != umig_cells.end(); ++cells_it)
				{
					int cell_with_same_umig_id = cells_it->first;
					if (cb_reassigned[cell_with_same_umig_id] == cur_gene.index)
						continue;

					if (this->_cells_genes[cb_reassigned[cell_with_same_umig_id]].size() > cur_gene.count)
					{
						umig_top[cb_reassigned[cell_with_same_umig_id]]++;
					}
				}
				umigs_count++;
			}
		}
		return umigs_count;
	}
}