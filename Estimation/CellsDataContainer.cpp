#include "CellsDataContainer.h"

#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"
#include "Tools/UtilFunctions.h"

#include "api/BamMultiReader.h"

#include <boost/range/adaptor/reversed.hpp>

using namespace std;

namespace Estimation
{
	const double EPS = 0.00001;

	CellsDataContainer::CellsDataContainer(double min_merge_fraction, int min_genes_before_merge, bool merge_tags,
	                                       int min_genes_after_merge, int max_merge_edit_distance, size_t top_print_size)
		: _merge_tags(merge_tags)
		, _top_print_size(top_print_size)
		, _min_merge_fraction(min_merge_fraction)
		, _min_genes_after_merge(min_genes_after_merge)
		, _min_genes_before_merge(min_genes_before_merge)
		, _max_merge_edit_distance(max_merge_edit_distance)
	{}

	void CellsDataContainer::init_filled(const CellsDataContainer::s_ii_hash_t & umig_cells_counts)
	{
		this->_cells_genes_counts_sorted = this->count_cells_genes();

		if (this->_merge_tags)
		{
			this->merge_genes(umig_cells_counts);

			this->_cells_genes_counts_sorted = this->count_cells_genes(false);
		}
		else
		{
			for (auto const &gene_count : boost::adaptors::reverse(this->_cells_genes_counts_sorted))
			{
				if (gene_count.count < this->_min_genes_after_merge)
					break;

				this->_filtered_cells.push_back(gene_count.index);
			}
		}
	}

	void CellsDataContainer::merge_genes(const s_ii_hash_t &umig_cells_counts)
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
		for (auto const &cur_gene : this->_cells_genes_counts_sorted)
		{ // iterate through the minimally-selected CBs, from low to high counts
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
			for (auto const &cell: umig_top)
			{
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

	bool CellsDataContainer::merge(int top_cell_ind, double top_cell_fraction, const IndexedCount &gene_count,
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
		for (auto const &gene: this->_cells_genes[cell_id])
		{
			for (auto const &umi_count: gene.second)
			{
				this->_cells_genes[top_cell_ind][gene.first][umi_count.first] += umi_count.second;
			}
		}

		this->reassign(cell_id, top_cell_ind, cb_reassigned, cb_reassigned_to);

		return true;
	}

	void CellsDataContainer::reassign(size_t cell_id, int target_cell_id, ints_t &cb_reassigned,
	                                  ISIHM &cb_reassigned_to) const
	{
		cb_reassigned[cell_id] = target_cell_id; // set reassignment mapping
		cb_reassigned_to[target_cell_id].insert(cell_id); // reassign current cell

		// transfer mapping of the cbs previously mapped to kid
		ISIHM::iterator k = cb_reassigned_to.find(cell_id);
		if (k == cb_reassigned_to.end())
			return;

		for (auto m: k->second)
		{
			cb_reassigned[m] = target_cell_id; // update reassignment mapping
			cb_reassigned_to[target_cell_id].insert(m);
		}

		k->second.clear();
		cb_reassigned_to[target_cell_id].insert(k->first); // reassign to the new cell
	}

	int CellsDataContainer::add_record(const string &cell_barcode, const string &umi, const string &gene, s_i_hash_t &cells_ids)
	{
		pair<s_i_hash_t::iterator, bool> res = cells_ids.emplace(cell_barcode, this->_cells_genes.size());
		if (res.second)
		{ // new cb
			this->_cells_genes.push_back(genes_t());
			this->_cells_names.push_back(cell_barcode);
		}
		int cell_id = res.first->second;
		this->_cells_genes[cell_id][gene][umi]++;

		L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi] << " gene=" <<
				this->_cells_genes[cell_id][gene].size()
				<< " CB=" << this->_cells_genes[cell_id].size();

		return cell_id;
	}

	CellsDataContainer::i_counter_t CellsDataContainer::count_cells_genes(bool logs) const
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

	string CellsDataContainer::get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const
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

	Stats& CellsDataContainer::stats()
	{
		return this->_stats;
	}

	const Stats& CellsDataContainer::stats() const
	{
		return this->_stats;
	}

	const CellsDataContainer::genes_t &CellsDataContainer::cell_genes(size_t index) const
	{
		return this->_cells_genes[index];
	}

	const CellsDataContainer::i_counter_t &CellsDataContainer::cells_genes_counts_sorted() const
	{
		return this->_cells_genes_counts_sorted;
	}

	const string &CellsDataContainer::cell_name(size_t index) const
	{
		return this->_cells_names[index];
	}

	const CellsDataContainer::ids_t &CellsDataContainer::filtered_cells() const
	{
		return this->_filtered_cells;
	}

	size_t CellsDataContainer::get_umig_top(const ints_t &cb_reassigned, const IndexedCount &cur_gene,
	                                        const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const
	{
		size_t umigs_count = 0;
		for (auto const &gene: this->_cells_genes[cur_gene.index])
		{
			const string &gene_name = gene.first;
			const s_i_hash_t &umis = gene.second;

			for (auto const &umi_count: umis)
			{
				string umig = umi_count.first + gene_name;
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
