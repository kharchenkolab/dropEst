#include "CellsDataContainer.h"

#include "Estimation/Merge/MergeStrategyAbstract.h"
#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"

#include "api/BamMultiReader.h"

using namespace std;

namespace Estimation
{
	using Tools::IndexedValue;

	CellsDataContainer::CellsDataContainer(std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy,
										   size_t top_print_size)
		: _merge_strategy(merge_strategy)
		, _top_print_size(top_print_size)
		, _is_initialized(false)
	{
		L_TRACE << this->_merge_strategy->merge_type() << " merge selected";
	}

	void CellsDataContainer::merge_and_filter()
	{
		if (!this->_is_initialized)
			throw runtime_error("You must initialize container");

		this->_merge_targets = this->_merge_strategy->merge(*this);
		this->stats().merge(this->_merge_targets, this->cell_barcodes_raw());
		
		this->update_cell_sizes(this->_merge_strategy->min_genes_after_merge());

		this->_filtered_cells.clear();
		for (auto const &val : this->cells_genes_counts_sorted())
		{
			this->_filtered_cells.push_back(val.index);
		}
	}

	size_t CellsDataContainer::add_record(const string &cell_barcode, const string &umi, const string &gene, s_i_map_t &cells_ids)
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
		size_t cell_id = res.first->second;
		this->_cells_genes[cell_id][gene][umi]++;

		L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi] << " gene=" <<
				this->_cells_genes[cell_id][gene].size() << " CB=" << this->_cells_genes[cell_id].size();

		return cell_id;
	}

	void CellsDataContainer::merge(size_t source_cell_ind, size_t target_cell_ind)
	{
		auto &target_cell = this->_cells_genes.at(target_cell_ind);
		for (auto const &gene: this->_cells_genes[source_cell_ind])
		{
			for (auto const &umi_count: gene.second)
			{
				target_cell[gene.first][umi_count.first] += umi_count.second;
			}
		}

		this->_is_cell_merged.at(source_cell_ind) = true;
	}

	void CellsDataContainer::exclude_cell(size_t index)
	{
		this->_is_cell_excluded.at(index) = true;
	}

	void CellsDataContainer::update_cell_sizes(int genes_threshold, bool logs)
	{
		this->_filtered_cells_genes_counts_sorted.clear(); // <genes_count,cell_id> pairs
		for (size_t i = 0; i < this->_cells_genes.size(); i++)
		{
			if (this->_is_cell_excluded[i] || this->_is_cell_merged[i])
				continue;

			size_t genes_count = this->_cells_genes[i].size();
			if (genes_count >= genes_threshold)
			{
				this->_filtered_cells_genes_counts_sorted.push_back(IndexedValue(i, genes_count));
			}
		}

		this->_cell_sizes.clear();
		this->_cell_sizes.resize(this->_cells_genes.size());
		for (size_t i = 0; i < this->_cells_genes.size(); ++i)
		{
			size_t cell_size = 0;
			for (auto const &gene: this->_cells_genes[i])
			{
				cell_size += gene.second.size();
			}
			this->_cell_sizes[i] = cell_size;
		}

		if (logs)
		{
			L_TRACE << this->_filtered_cells_genes_counts_sorted.size() << " CBs with more than " << genes_threshold << " genes";
		}

		sort(this->_filtered_cells_genes_counts_sorted.begin(), this->_filtered_cells_genes_counts_sorted.end(), IndexedValue::value_less);

		if (logs)
		{
			L_TRACE << this->get_cb_count_top_verbose(this->_filtered_cells_genes_counts_sorted);
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
		return this->_cells_genes.at(index);
	}

	const CellsDataContainer::i_counter_t &CellsDataContainer::cells_genes_counts_sorted() const
	{
		return this->_filtered_cells_genes_counts_sorted;
	}

	const string &CellsDataContainer::cell_barcode(size_t index) const
	{
		return this->_cell_barcodes.at(index);
	}

	const CellsDataContainer::names_t& CellsDataContainer::cell_barcodes_raw() const
	{
		return this->_cell_barcodes;
	}

	const CellsDataContainer::ids_t &CellsDataContainer::filtered_cells() const
	{
		return this->_filtered_cells;
	}

	void CellsDataContainer::set_initialized()
	{
		if (this->_is_initialized)
			throw std::runtime_error("Container is already initialized");

		this->_is_cell_excluded = flags_t(this->_cell_barcodes.size(), false);
		this->_is_cell_merged = flags_t(this->_cell_barcodes.size(), false);

		this->update_cell_sizes(this->_merge_strategy->min_genes_before_merge());
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

	CellsDataContainer::s_ul_hash_t CellsDataContainer::umis_distribution() const
	{
		s_ul_hash_t umis_dist;
		for (auto const &cell : this->_filtered_cells_genes_counts_sorted)
		{
			for (auto const &gene : this->_cells_genes[cell.index])
			{
				for (auto const &umi : gene.second)
				{
					umis_dist[umi.first]++;
				}
			}
		}

		return umis_dist;
	}

	bool CellsDataContainer::is_cell_merged(size_t cell_id) const
	{
		return this->_is_cell_merged.at(cell_id);
	}

	std::string CellsDataContainer::merge_type() const
	{
		return this->_merge_strategy->merge_type();
	}

	const size_t CellsDataContainer::cell_size(size_t cell_index) const
	{
		return this->_cell_sizes.at(cell_index);
	}

	const CellsDataContainer::ids_t &CellsDataContainer::merge_targets() const
	{
		return this->_merge_targets;
	}

	bool CellsDataContainer::is_cell_excluded(size_t cell_id) const
	{
		return this->_is_cell_excluded.at(cell_id);
	}
}
