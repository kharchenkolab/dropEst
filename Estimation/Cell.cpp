#include "Cell.h"

#include <cstring>

namespace Estimation
{
	Cell::Cell(const std::string &barcode, size_t min_genes_to_be_real, StringIndexer *gene_indexer,
	           StringIndexer *umi_indexer)
		: _barcode(std::unique_ptr<char[]>(new char[barcode.length() + 1], std::default_delete<char[]>()))
		, _min_genes_to_be_real(min_genes_to_be_real)
		, _is_merged(false)
		, _is_excluded(false)
		, _requested_genes_num(0)
		, _requested_umis_num(0)
		, _gene_indexer(gene_indexer)
		, _umi_indexer(umi_indexer)
	{
		strcpy(this->_barcode.get(), barcode.c_str());
	}

	void Cell::add_umi(const std::string &gene, const std::string &umi, const UMI::Mark &umi_mark)
	{
		auto gene_it = this->_genes.emplace(this->_gene_indexer->add(gene), this->_umi_indexer);
		bool is_new = gene_it.first->second.add_umi(umi, umi_mark);
		if (is_new)
		{
			this->_stats.inc(Stats::TOTAL_UMIS_PER_CB);
		}
	}

	void Cell::set_merged()
	{
		this->_is_merged = true;
	}

	void Cell::set_excluded()
	{
		this->_is_excluded = true;
	}

	void Cell::merge(const Cell &source)
	{
		for (auto const &gene: source._genes)
		{
			auto gene_it = this->_genes.emplace(gene.first, this->_umi_indexer);
			gene_it.first->second.merge(gene.second);
		}

		this->_stats.merge(source.stats());
	}

	void Cell::merge_umis(StringIndexer::index_t gene, const s_s_hash_t &merge_targets)
	{
		auto &gene_umis = this->_genes.at(gene);
		for (auto const &target: merge_targets)
		{
			if (target.second == target.first)
				continue;

			gene_umis.merge(target.first, target.second);
			this->_stats.dec(Stats::TOTAL_UMIS_PER_CB);
		}
	}

	const Cell::genes_t &Cell::genes() const
	{
		return this->_genes;
	}

	Cell::s_ul_hash_t Cell::requested_umis_per_gene(const UMI::Mark::query_t &query_marks, bool return_reads) const
	{
		s_ul_hash_t umis_per_gene;
		for (auto const &gene : this->_genes)
		{
			size_t umis_num = gene.second.number_of_requested_umis(query_marks, return_reads);

			if (umis_num == 0)
				continue;

			umis_per_gene.emplace(this->_gene_indexer->get_value(gene.first), umis_num);
		}

		return umis_per_gene;
	}

	Cell::ss_ul_hash_t Cell::requested_reads_per_umi_per_gene(const UMI::Mark::query_t &query_marks) const
	{
		ss_ul_hash_t reads_per_umi_per_gene;
		for (auto const &gene : this->_genes)
		{
			s_ul_hash_t reads_per_umi(gene.second.requested_reads_per_umi(query_marks));
			if (reads_per_umi.empty())
				continue;

			reads_per_umi_per_gene.emplace(this->_gene_indexer->get_value(gene.first), reads_per_umi);
		}

		return reads_per_umi_per_gene;
	}

	bool Cell::is_merged() const
	{
		return this->_is_merged;
	}

	bool Cell::is_excluded() const
	{
		return this->_is_excluded;
	}

	std::string Cell::barcode() const
	{
		return std::string(this->_barcode.get());
	}

	const char* Cell::barcode_c() const
	{
		return this->_barcode.get();
	}

	size_t Cell::umis_number() const
	{
		return this->stats().get(Stats::TOTAL_UMIS_PER_CB);
	}

	size_t Cell::requested_genes_num() const
	{
		return this->_requested_genes_num;
	}

	size_t Cell::requested_umis_num() const
	{
		return this->_requested_umis_num;
	}

	size_t Cell::size() const
	{
		return this->_genes.size();
	}

	bool Cell::is_real() const
	{
		return !this->_is_excluded && !this->_is_merged && this->size() >= this->_min_genes_to_be_real;
	}

	void Cell::update_requested_size(const UMI::Mark::query_t &query_marks)
	{
		this->_requested_genes_num = 0;
		this->_requested_umis_num = 0;
		for (auto const &gene : this->_genes)
		{
			size_t cur_num = gene.second.number_of_requested_umis(query_marks, false);
			if (cur_num == 0)
				continue;

			this->_requested_umis_num += cur_num;
			this->_requested_genes_num++;
		}
	}

	const Gene& Cell::at(const std::string &gene) const
	{
		return this->_genes.at(this->_gene_indexer->get_index(gene));
	}

	const Stats &Cell::stats() const
	{
		return this->_stats;
	}

	Stats &Cell::stats()
	{
		return this->_stats;
	}
}
