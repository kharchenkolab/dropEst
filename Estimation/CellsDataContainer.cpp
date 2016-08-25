#include "CellsDataContainer.h"

#include "Estimation/Merge/IMergeStrategy.h"
#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"

#include "api/BamMultiReader.h"

using namespace std;

namespace Estimation
{
	using Tools::IndexedValue;

	CellsDataContainer::CellsDataContainer(std::shared_ptr<Merge::IMergeStrategy> merge_strategy,
										   size_t top_print_size)
		: _merge_strategy(merge_strategy)
		, _top_print_size(top_print_size)
		, _is_initialized(false)
	{}

	void CellsDataContainer::merge_and_filter(const s_uu_hash_t &umig_cells_counts)
	{
		if (!this->_is_initialized)
			throw runtime_error("You must initialize container");

		this->_merge_strategy->merge(*this, umig_cells_counts, this->_filtered_cells);
		this->update_cells_genes_counts(this->_merge_strategy->min_genes_after_merge(), false);
	}

	int CellsDataContainer::add_record(const string &cell_barcode, const string &umi, const string &gene, s_i_map_t &cells_ids)
	{
		if (this->_is_initialized)
			throw runtime_error("Container is already initialized");

		pair<s_i_map_t::iterator, bool> res = cells_ids.emplace(cell_barcode, this->_cells_genes.size());
		if (res.second)
		{ // new cb
			this->_cells_genes.push_back(genes_t());
			this->_cell_ids_by_cb[cell_barcode] = this->_cell_barcodes.size();
			this->_cell_barcodes.push_back(cell_barcode);
		}
		int cell_id = res.first->second;
		this->_cells_genes[cell_id][gene][umi]++;

		L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi] << " gene=" <<
				this->_cells_genes[cell_id][gene].size() << " CB=" << this->_cells_genes[cell_id].size();

		return cell_id;
	}

	void CellsDataContainer::merge(size_t source_cell_ind, size_t target_cell_ind)
	{
		for (auto const &gene: this->_cells_genes[source_cell_ind])
		{
			for (auto const &umi_count: gene.second)
			{
				this->_cells_genes[target_cell_ind][gene.first][umi_count.first] += umi_count.second;
			}
		}
	}

	void CellsDataContainer::exclude_cell(size_t index)
	{
		this->_is_cell_excluded[index] = true;
	}

	void CellsDataContainer::update_cells_genes_counts(int threshold, bool logs)
	{
		this->filtered_cells_genes_counts_sorted.clear(); // <genes_count,cell_id> pairs
		for (size_t i = 0; i < this->_cells_genes.size(); i++)
		{
			if (this->_is_cell_excluded[i])
				continue;

			size_t genes_count = this->_cells_genes[i].size();
			if (genes_count >= threshold)
			{
				this->filtered_cells_genes_counts_sorted.push_back(IndexedValue(i, genes_count));
			}
		}

		if (logs)
		{
			L_TRACE << this->filtered_cells_genes_counts_sorted.size() << " CBs with more than " << threshold << " genes";
		}

		sort(this->filtered_cells_genes_counts_sorted.begin(), this->filtered_cells_genes_counts_sorted.end(), IndexedValue::value_less);

		if (logs)
		{
			L_TRACE << this->get_cb_count_top_verbose(this->filtered_cells_genes_counts_sorted);
		}
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
				ss << cells_genes_counts[i].value << "\t" << this->_cell_barcodes[cells_genes_counts[i].index] << "\n";
			}
		}
		else
		{
			ss << "no valid CBs found\n";
		}

		return ss.str();
	}

	Stats& CellsDataContainer::stats() const
	{
		return this->_stats;
	}

	const CellsDataContainer::genes_t &CellsDataContainer::cell_genes(size_t index) const
	{
		return this->_cells_genes[index];
	}

	const CellsDataContainer::i_counter_t &CellsDataContainer::cells_genes_counts_sorted() const
	{
		return this->filtered_cells_genes_counts_sorted;
	}

	const string &CellsDataContainer::cell_barcode(size_t index) const
	{
		return this->_cell_barcodes[index];
	}

	const CellsDataContainer::names_t& CellsDataContainer::cell_barcodes() const
	{
		return this->_cell_barcodes;
	}

	const CellsDataContainer::ids_t &CellsDataContainer::filtered_cells() const
	{
		return this->_filtered_cells;
	}

	void CellsDataContainer::set_initialized()
	{
		this->_is_cell_excluded = flags_t(this->_cell_barcodes.size(), false);
		this->update_cells_genes_counts(this->_merge_strategy->min_genes_before_merge());
		this->_is_initialized = true;
	}

	CellsDataContainer::names_t CellsDataContainer::excluded_cells() const
	{
		names_t ec;
		for (size_t i = 0; i < this->_is_cell_excluded.size(); ++i)
		{
			if (this->_is_cell_excluded[i])
			{
				ec.push_back(this->_cell_barcodes[i]);
			}
		}

		return ec;
	}

	const CellsDataContainer::s_ul_hash_t& CellsDataContainer::cell_ids_by_cb() const
	{
		return this->_cell_ids_by_cb;
	}
}
