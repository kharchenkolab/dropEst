#include "CellsDataContainer.h"

#include <Estimation/Merge/MergeStrategyAbstract.h>
#include <Estimation/Merge/UMIs/MergeUMIsStrategySimple.h>

#include <Tools/Logs.h>
#include <Tools/RefGenesContainer.h>

#include <api/BamMultiReader.h>

using namespace std;

namespace Estimation
{
	using Tools::IndexedValue;

	const std::string UMI::Mark::DEFAULT_CODE = "eEBA";
	const size_t CellsDataContainer::TOP_PRINT_SIZE = 10;

	CellsDataContainer::CellsDataContainer(std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy,
	                                       std::shared_ptr<Merge::UMIs::MergeUMIsStrategySimple> umi_merge_strategy,
		                                   const std::vector<UMI::Mark> &gene_match_levels, int max_cells_num)
		: _merge_strategy(merge_strategy)
		, _umi_merge_strategy(umi_merge_strategy)
		, _max_cells_num(max_cells_num)
		, _is_initialized(false)
		, _query_marks(gene_match_levels)
		, _has_exon_reads(0)
		, _has_intron_reads(0)
		, _has_not_annotated_reads(0)
		, _number_of_real_cells(0)
	{
		L_TRACE << this->_merge_strategy->merge_type() << " merge selected.";
	}

	void CellsDataContainer::merge_and_filter()
	{
		if (!this->_is_initialized)
			throw runtime_error("You must initialize container");

		this->_merge_targets = this->_merge_strategy->merge(*this);
		this->stats().merge(this->_merge_targets, this->_cells);

		this->_umi_merge_strategy->merge(*this);
		size_t filtered_cells_num = this->update_cell_sizes(this->_merge_strategy->min_genes_after_merge(), this->_max_cells_num);

		L_TRACE << this->_number_of_real_cells << " cells are considered as real.\n";

		L_TRACE << filtered_cells_num << " CBs with more than " << this->_merge_strategy->min_genes_after_merge()
		        << " genes, which have UMIs of the requested type.";

		L_TRACE << this->get_cb_count_top_verbose(this->_filtered_gene_counts_sorted);

		this->_filtered_cells.clear();
		for (auto const &val : this->cells_gene_counts_sorted())
		{
			this->_filtered_cells.push_back(val.index);
		}
	}

	size_t CellsDataContainer::add_record(const string &cell_barcode, const string &umi, const string &gene,
	                                      const UMI::Mark &umi_mark)
	{
		if (this->_is_initialized)
			throw runtime_error("Container is already initialized");

		if (umi_mark.check(UMI::Mark::HAS_EXONS))
		{
			++this->_has_exon_reads;
		}
		if (umi_mark.check(UMI::Mark::HAS_INTRONS))
		{
			++this->_has_intron_reads;
		}
		if (umi_mark.check(UMI::Mark::HAS_NOT_ANNOTATED))
		{
			++this->_has_not_annotated_reads;
		}

		if (umi_mark == UMI::Mark::NONE)
		{
			L_WARN << "Empty mark for CB '" + cell_barcode + "', UMI '" + umi + "', gene '" + gene + "'";
		}

		auto res = this->_cell_ids_by_cb.emplace(cell_barcode, this->_cell_ids_by_cb.size());
		if (res.second)
		{ // new cb
			this->_cells.push_back(Cell(cell_barcode, this->_merge_strategy->min_genes_before_merge(), this->_query_marks));
		}
		size_t cell_id = res.first->second;
		this->_cells[cell_id].add_umi(gene, umi, umi_mark);

		return cell_id;
	}

	void CellsDataContainer::merge_cells(size_t source_cell_ind, size_t target_cell_ind)
	{
		this->_cells.at(target_cell_ind).merge(this->_cells.at(source_cell_ind));
		this->_cells.at(source_cell_ind).set_merged();
	}

	void CellsDataContainer::exclude_cell(size_t index)
	{
		this->_cells.at(index).set_excluded();
	}

	size_t CellsDataContainer::update_cell_sizes(size_t requested_genes_threshold, int cell_threshold)
	{
		this->_number_of_real_cells = 0;
		for (size_t i = 0; i < this->_cells.size(); i++)
		{
			this->_cells[i].update_requested_size();
			if (!this->_cells[i].is_real())
				continue;

			this->_number_of_real_cells++;
		}

		return this->get_filtered_gene_counts(requested_genes_threshold, cell_threshold, this->_filtered_gene_counts_sorted);
	}

	string CellsDataContainer::get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const
	{
		stringstream ss;
		if (cells_genes_counts.size() > 0)
		{
			ss << "top CBs:\n";
			long low_border = cells_genes_counts.size() - min(cells_genes_counts.size(), CellsDataContainer::TOP_PRINT_SIZE - 1);
			for (long i = cells_genes_counts.size() - 1; i >= low_border; --i)
			{
				ss << cells_genes_counts[i].value << "\t" << this->_cells[cells_genes_counts[i].index].barcode() << "\n";
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

	const Cell &CellsDataContainer::cell(size_t index) const
	{
		return this->_cells.at(index);
	}

	const CellsDataContainer::i_counter_t &CellsDataContainer::cells_gene_counts_sorted() const
	{
		return this->_filtered_gene_counts_sorted;
	}

	const CellsDataContainer::ids_t &CellsDataContainer::filtered_cells() const
	{
		return this->_filtered_cells;
	}

	void CellsDataContainer::set_initialized()
	{
		if (this->_is_initialized)
			throw std::runtime_error("Container is already initialized");

		this->update_cell_sizes(0, -1);

		L_TRACE << "\n" << this->_filtered_gene_counts_sorted.size() << " CBs with more than "
		        << this->_merge_strategy->min_genes_before_merge() << " genes";
		L_TRACE << this->get_cb_count_top_verbose(this->_filtered_gene_counts_sorted);

		this->_is_initialized = true;
	}

	CellsDataContainer::names_t CellsDataContainer::excluded_cells() const
	{
		names_t ec;
		for (auto const &cell : this->_cells)
		{
			if (cell.is_excluded())
			{
				ec.push_back(cell.barcode());
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
		for (auto const &cell : this->_filtered_gene_counts_sorted)
		{
			for (auto const &gene : this->_cells[cell.index].genes())
			{
				for (auto const &umi : gene.second.umis())
				{
					umis_dist[umi.first]++;
				}
			}
		}

		return umis_dist;
	}

	std::string CellsDataContainer::merge_type() const
	{
		return this->_merge_strategy->merge_type();
	}

	const CellsDataContainer::ids_t &CellsDataContainer::merge_targets() const
	{
		return this->_merge_targets;
	}

	void CellsDataContainer::merge_umis(size_t cell_id, const std::string &gene,
	                                    const CellsDataContainer::s_s_hash_t &merge_targets)
	{
		this->_cells.at(cell_id).merge_umis(gene, merge_targets);
	}

	const std::vector<UMI::Mark>& CellsDataContainer::gene_match_level() const
	{
		return this->_query_marks;
	}

	size_t CellsDataContainer::has_exon_reads_num() const
	{
		return this->_has_exon_reads;
	}

	size_t CellsDataContainer::has_intron_reads_num() const
	{
		return this->_has_intron_reads;
	}

	size_t CellsDataContainer::has_not_annotated_reads_num() const
	{
		return this->_has_not_annotated_reads;
	}

	size_t CellsDataContainer::real_cells_number() const
	{
		return this->_number_of_real_cells;
	}

	size_t CellsDataContainer::total_cells_number() const
	{
		return this->_cells.size();
	}

	size_t CellsDataContainer::get_filtered_gene_counts(size_t requested_genes_threshold, int cell_threshold,
	                                                    i_counter_t &filtered_gene_counts_sorted) const
	{
		filtered_gene_counts_sorted.clear();
		for (size_t i = 0; i < this->_cells.size(); i++)
		{
			if (!this->_cells[i].is_real())
				continue;

			size_t genes_count = this->_cells[i].requested_genes_num();
			if (genes_count >= requested_genes_threshold)
			{
				filtered_gene_counts_sorted.push_back(IndexedValue(i, genes_count));
			}
		}

		sort(filtered_gene_counts_sorted.begin(), filtered_gene_counts_sorted.end(), IndexedValue::value_less);

		size_t filtered_cells_num = filtered_gene_counts_sorted.size();

		if (cell_threshold > 0 && cell_threshold < filtered_gene_counts_sorted.size())
		{
			filtered_gene_counts_sorted.erase(filtered_gene_counts_sorted.begin(),
			                                  filtered_gene_counts_sorted.end() - unsigned(cell_threshold));
		}

		return filtered_cells_num;
	}
}
