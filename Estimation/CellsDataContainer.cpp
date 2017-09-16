#include "CellsDataContainer.h"

#include <Estimation/Merge/MergeStrategyAbstract.h>
#include <Estimation/Merge/UMIs/MergeUMIsStrategySimple.h>

#include <Tools/Logs.h>
#include <Tools/RefGenesContainer.h>

#include <api/BamMultiReader.h>

#include <numeric>
#include <boost/bind.hpp>

using namespace std;

namespace Estimation
{
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

		this->_umi_merge_strategy->merge(*this);
		size_t filtered_cells_num = this->update_cell_sizes(this->_merge_strategy->min_genes_after_merge(), this->_max_cells_num);

		L_TRACE << this->_number_of_real_cells << " cells are considered as real.\n";

		L_TRACE << filtered_cells_num << " CBs with more than " << this->_merge_strategy->min_genes_after_merge()
		        << " genes, which have UMIs of the requested type.";

		L_TRACE << this->get_cb_count_top_verbose();
	}

	void CellsDataContainer::add_record(const ReadInfo &read_info)
	{
		if (this->_is_initialized)
			throw runtime_error("Container is already initialized");

		auto res = this->_cell_ids_by_cb.emplace(read_info.params.cell_barcode(), this->_cell_ids_by_cb.size());
		if (res.second)
		{
			this->_cells.push_back(Cell(read_info.params.cell_barcode(), this->_merge_strategy->min_genes_before_merge(),
			                            this->_query_marks, &this->_gene_indexer, &this->_umi_indexer));
		}

		size_t cell_id = res.first->second;

		if (read_info.gene.empty())
		{
			this->_cells[cell_id].stats().inc(Stats::INTERGENIC_READS_PER_CHR_PER_CELL, read_info.chromosome_name);
			return;
		}

		if (read_info.umi_mark == UMI::Mark::NONE)
		{
			L_WARN << "Empty mark for CB '" + read_info.params.cell_barcode() + "', UMI '" + read_info.params.umi() +
						"', gene '" + read_info.gene + "'";
		}

		this->_cells[cell_id].add_umi(read_info);
		this->update_cell_stats(cell_id, read_info.umi_mark, read_info.chromosome_name);
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

		return this->update_filtered_gene_counts(requested_genes_threshold, cell_threshold);
	}

	string CellsDataContainer::get_cb_count_top_verbose() const
	{
		stringstream ss;
		if (this->_filtered_cells.size() > 0)
		{
			ss << "top CBs:\n";
			long low_border = this->_filtered_cells.size() - min(this->_filtered_cells.size(), CellsDataContainer::TOP_PRINT_SIZE - 1);
			for (long i = this->_filtered_cells.size() - 1; i >= low_border; --i)
			{
				auto const &cur_cell = this->_cells[this->_filtered_cells[i]];
				ss << cur_cell.requested_genes_num() << "\t" << cur_cell.barcode() << "\n";
			}
		}
		else
		{
			ss << "no valid CBs found\n";
		}

		return ss.str();
	}

	const Cell &CellsDataContainer::cell(size_t index) const
	{
		return this->_cells.at(index);
	}

	Cell &CellsDataContainer::cell(size_t index)
	{
		return this->_cells.at(index);
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

		L_TRACE << "\n" << this->_filtered_cells.size() << " CBs with more than "
		        << this->_merge_strategy->min_genes_before_merge() << " genes";
		L_TRACE << this->get_cb_count_top_verbose();

		this->_is_initialized = true;
	}

	size_t CellsDataContainer::cell_id_by_cb(const std::string &barcode) const
	{
		return this->_cell_ids_by_cb.at(barcode);
	}

	CellsDataContainer::s_ul_hash_t CellsDataContainer::umis_distribution() const
	{
		s_ul_hash_t umis_dist;
		for (size_t cell_id : this->_filtered_cells)
		{
			for (auto const &gene : this->_cells[cell_id].genes())
			{
				for (auto const &umi : gene.second.umis())
				{
					umis_dist[this->_umi_indexer.get_value(umi.first)]++;
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

	void CellsDataContainer::merge_umis(size_t cell_id, StringIndexer::index_t gene,
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

	size_t CellsDataContainer::update_filtered_gene_counts(size_t requested_genes_threshold, int cell_threshold)
	{
		this->_filtered_cells.clear();
		for (size_t i = 0; i < this->_cells.size(); i++)
		{
			if (!this->_cells[i].is_real())
				continue;

			size_t genes_count = this->_cells[i].requested_genes_num();
			if (genes_count >= requested_genes_threshold)
			{
				this->_filtered_cells.push_back(i);
			}
		}

		sort(this->_filtered_cells.begin(), this->_filtered_cells.end(), boost::bind(&CellsDataContainer::compare_cells, this, _1, _2));

		size_t filtered_cells_num = this->_filtered_cells.size();

		if (cell_threshold > 0 && cell_threshold < this->_filtered_cells.size())
		{
			this->_filtered_cells.erase(this->_filtered_cells.begin(),
			                            this->_filtered_cells.end() - unsigned(cell_threshold));
		}

		return filtered_cells_num;
	}

	CellsDataContainer::s_i_hash_t CellsDataContainer::get_stat_by_real_cells(Stats::CellStatType type) const
	{
		s_i_hash_t res;
		for (auto const &cell : this->_cells)
		{
			if (!cell.is_real())
				continue;

			res[cell.barcode()] = cell.stats().get(type);
		}

		return res;
	}

	void CellsDataContainer::get_stat_by_real_cells(Stats::CellChrStatType stat, names_t &cell_barcodes,
	                                                names_t &chromosome_names, counts_t &counts) const
	{
		for (auto const &cell : this->_cells)
		{
			if (!cell.is_real())
				continue;

			if (cell.stats().get(stat, counts))
			{
				cell_barcodes.push_back(cell.barcode());
			}
		}

		chromosome_names = Stats::presented_chromosomes(stat);
	}

	void CellsDataContainer::update_cell_stats(size_t cell_id, const UMI::Mark &mark, const std::string &chromosome_name)
	{
		auto &cell = this->_cells[cell_id];
		cell.stats().inc(Stats::TOTAL_READS_PER_CB);
		if (mark.check(UMI::Mark::HAS_EXONS))
		{
			cell.stats().inc(Stats::EXON_READS_PER_CHR_PER_CELL, chromosome_name);
			++this->_has_exon_reads;
		}
		if (mark.check(UMI::Mark::HAS_INTRONS))
		{
			cell.stats().inc(Stats::INTRON_READS_PER_CHR_PER_CELL, chromosome_name);
			++this->_has_intron_reads;
		}
		if (mark.check(UMI::Mark::HAS_NOT_ANNOTATED))
		{
			++this->_has_not_annotated_reads;
		}
	}

	bool CellsDataContainer::compare_cells(size_t cell1_id, size_t cell2_id) const
	{
		auto const &cell1 = this->_cells[cell1_id];
		auto const &cell2 = this->_cells[cell2_id];

		if (cell1.requested_genes_num() != cell2.requested_genes_num())
			return cell1.requested_genes_num() < cell2.requested_genes_num();

		if (cell1.requested_umis_num() != cell2.requested_umis_num())
			return cell1.requested_umis_num() < cell2.requested_umis_num();

		if (cell1.umis_number() != cell2.umis_number())
			return cell1.umis_number() < cell2.umis_number();

		return cell1.barcode() < cell2.barcode();
	}

	const StringIndexer &CellsDataContainer::gene_indexer() const
	{
		return this->_gene_indexer;
	}

	const StringIndexer &CellsDataContainer::umi_indexer() const
	{
		return this->_umi_indexer;
	}
}
