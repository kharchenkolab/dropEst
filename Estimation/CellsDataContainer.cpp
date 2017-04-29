#include "CellsDataContainer.h"

#include <Estimation/Merge/MergeStrategyAbstract.h>
#include <Estimation/MergeUMIs/MergeUMIsStrategySimple.h>

#include <Tools/Logs.h>
#include <Tools/RefGenesContainer.h>

#include <api/BamMultiReader.h>

using namespace std;

namespace Estimation
{
	using Tools::IndexedValue;

	CellsDataContainer::CellsDataContainer(std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy,
	                                       std::shared_ptr<MergeUMIs::MergeUMIsStrategySimple> umi_merge_strategy,
	                                       size_t top_print_size, GeneMatchLevel gene_match_level)
		: _merge_strategy(merge_strategy)
		, _umi_merge_strategy(umi_merge_strategy)
		, _top_print_size(top_print_size)
		, _is_initialized(false)
		, _gene_match_level(gene_match_level)
	{
		L_TRACE << this->_merge_strategy->merge_type() << " merge selected";
	}

	void CellsDataContainer::merge_and_filter()
	{
		if (!this->_is_initialized)
			throw runtime_error("You must initialize container");

		this->_merge_targets = this->_merge_strategy->merge(*this);
		this->stats().merge(this->_merge_targets, this->cell_barcodes_raw());

//		this->_umi_merge_strategy->merge(*this); TODO: uncomment and debug
		this->remove_excluded_umis();
		this->update_cell_sizes(this->_merge_strategy->min_genes_after_merge());

		this->_filtered_cells.clear();
		for (auto const &val : this->cells_gene_counts_sorted())
		{
			this->_filtered_cells.push_back(val.index);
		}
	}

	size_t CellsDataContainer::add_record(const string &cell_barcode, const string &umi, const string &gene,
	                                      const Mark &umi_mark)
	{
		if (this->_is_initialized)
			throw runtime_error("Container is already initialized");

		auto res = this->_cell_ids_by_cb.emplace(cell_barcode, this->_cell_barcodes.size());
		if (res.second)
		{ // new cb
			this->_cells_genes.push_back(genes_t());
			this->_cell_barcodes.push_back(cell_barcode);
		}
		size_t cell_id = res.first->second;
		auto & new_umi = this->_cells_genes[cell_id][gene].emplace(umi, 0).first->second;
		new_umi.read_count++;
		new_umi.mark.add(umi_mark);

		L_DEBUG << "CB/UMI=" << this->_cells_genes[cell_id][gene][umi].read_count << " gene=" <<
		        this->_cells_genes[cell_id][gene].size() << " CB=" << this->_cells_genes[cell_id].size();

		return cell_id;
	}

	void CellsDataContainer::merge_cells(size_t source_cell_ind, size_t target_cell_ind)
	{
		auto &target_cell = this->_cells_genes.at(target_cell_ind);
		for (auto const &gene: this->_cells_genes.at(source_cell_ind))
		{
			auto &target_gene = target_cell[gene.first];
			for (auto const &merged_umi: gene.second)
			{
				target_gene[merged_umi.first].merge(merged_umi.second);
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
		this->_filtered_cells_gene_counts_sorted.clear(); // <genes_count,cell_id> pairs
		for (size_t i = 0; i < this->_cells_genes.size(); i++)
		{
			if (this->_is_cell_excluded[i] || this->_is_cell_merged[i])
				continue;

			size_t genes_count = this->_cells_genes[i].size();
			if (genes_count >= genes_threshold)
			{
				this->_filtered_cells_gene_counts_sorted.push_back(IndexedValue(i, genes_count));
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
			L_TRACE << this->_filtered_cells_gene_counts_sorted.size() << " CBs with more than " << genes_threshold << " genes";
		}

		sort(this->_filtered_cells_gene_counts_sorted.begin(), this->_filtered_cells_gene_counts_sorted.end(), IndexedValue::value_less);

		if (logs)
		{
			L_TRACE << this->get_cb_count_top_verbose(this->_filtered_cells_gene_counts_sorted);
		}
	}

	string CellsDataContainer::get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const
	{
		stringstream ss;
		if (cells_genes_counts.size() > 0)
		{
			ss << "top CBs:\n";
			long low_border = cells_genes_counts.size() - min(cells_genes_counts.size(), this->_top_print_size - 1);
			for (long i = cells_genes_counts.size() - 1; i >= low_border; --i)
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

	const CellsDataContainer::i_counter_t &CellsDataContainer::cells_gene_counts_sorted() const
	{
		return this->_filtered_cells_gene_counts_sorted;
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
		for (auto const &cell : this->_filtered_cells_gene_counts_sorted)
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

	size_t CellsDataContainer::cell_size(size_t cell_index) const
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

	void CellsDataContainer::remove_excluded_umis()
	{
		for (auto &cell  : this->_cells_genes)
		{
			auto gene_iter = cell.begin();
			while (gene_iter != cell.end())
			{
				auto umi_iter = gene_iter->second.begin();
				while (umi_iter != gene_iter->second.end())
				{
					if (!umi_iter->second.mark.match(this->_gene_match_level))
					{
						umi_iter = gene_iter->second.erase(umi_iter);
					}
					else
					{
						++umi_iter;
					}
				}

				if (gene_iter->second.empty())
				{
					gene_iter = cell.erase(gene_iter);
				}
				else
				{
					++gene_iter;
				}
			}
		}
	}

	void CellsDataContainer::merge_umis(size_t cell_id, const std::string &gene,
	                                    const CellsDataContainer::s_s_hash_t &merge_targets)
	{
		auto &gene_umis = this->_cells_genes.at(cell_id).at(gene);
		for (auto const &target: merge_targets)
		{
			if (target.second == target.first)
				continue;

			gene_umis[target.second].merge(gene_umis.at(target.first));
			gene_umis.erase(target.first);
		}
	}

	CellsDataContainer::GeneMatchLevel CellsDataContainer::gene_match_level() const
	{
		return this->_gene_match_level;
	}

	CellsDataContainer::Mark::Mark(Mark::MarkType type)
			: _mark(type)
	{}

	void CellsDataContainer::Mark::add(Mark::MarkType type)
	{
		this->_mark |= type;
	}

	bool CellsDataContainer::Mark::check(Mark::MarkType type) const
	{
		return this->_mark & type;
	}

	void CellsDataContainer::Mark::add(const CellsDataContainer::Mark &mark)
	{
		this->_mark |= mark._mark;
	}

	bool CellsDataContainer::Mark::match(GeneMatchLevel match_level) const
	{
		switch (match_level)
		{
			case ANY:
				return this->check(HAS_EXONS);
			case BOTH_EXON:
				return !this->check(HAS_NOT_ANNOTATED) & !this->check(HAS_INTRONS)  & this->check(HAS_EXONS);
			case INTRON_EXON:
				return this->check(HAS_INTRONS) & this->check(HAS_EXONS);
			default:
				throw std::runtime_error("Unexpected gene match level: " + std::to_string(match_level));
		}
	}

	void CellsDataContainer::Mark::add(Tools::GtfRecord::RecordType type)
	{
		switch (type)
		{
			case Tools::GtfRecord::EXON:
				this->add(HAS_EXONS);
				break;
			case Tools::GtfRecord::INTRON:
				this->add(HAS_INTRONS);
				break;
			default:
				throw std::runtime_error("Unexpected GtfRecord type: " + std::to_string(type));
		}
	}

	bool CellsDataContainer::Mark::operator==(const CellsDataContainer::Mark::MarkType &other) const
	{
		return this->_mark == other;
	}

	bool CellsDataContainer::Mark::operator==(const CellsDataContainer::Mark &other) const
	{
		return this->_mark == other._mark;
	}

	CellsDataContainer::UMI::UMI(size_t read_count)
		: read_count(read_count)
		, mark()
	{}

	void CellsDataContainer::UMI::merge(const CellsDataContainer::UMI &umi)
	{
		this->read_count += umi.read_count;
		this->mark.add(umi.mark);
	}
}
