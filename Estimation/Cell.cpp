#include "Cell.h"

namespace Estimation
{
	Cell::Cell(const std::string &barcode, size_t min_genes_to_be_real, const std::vector<UMI::Mark> &query_marks)
			: _barcode(barcode)
			, _min_genes_to_be_real(min_genes_to_be_real)
			, _query_marks(query_marks)
			, _is_merged(false)
			, _is_excluded(false)
			, _umis_number(0)
			, _requested_genes_num(0)
	{}

	void Cell::add_umi(const std::string &gene, const std::string &umi, const UMI::Mark &umi_mark)
	{
		auto insert_it = this->_genes[gene].emplace(umi, 0);
		auto & new_umi = insert_it.first->second;
		new_umi.read_count++;
		new_umi.mark.add(umi_mark);

		if (insert_it.second)
		{
			this->_umis_number++;
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
			auto &target_gene = this->_genes[gene.first];
			for (auto const &merged_umi: gene.second)
			{
				auto insert_it = target_gene.insert(merged_umi);
				if (insert_it.second)
				{
					this->_umis_number++;
					continue;
				}
				insert_it.first->second.merge(merged_umi.second);
			}
		}
	}

	void Cell::merge_umis(const std::string &gene, const s_s_hash_t &merge_targets)
	{
		auto &gene_umis = this->_genes.at(gene);
		for (auto const &target: merge_targets)
		{
			if (target.second == target.first)
				continue;

			gene_umis[target.second].merge(gene_umis.at(target.first));
			gene_umis.erase(target.first);
		}
	}

	const Cell::genes_t &Cell::genes() const
	{
		return this->_genes;
	}

	Cell::s_ul_hash_t Cell::requested_umis_per_gene(bool return_reads) const
	{
		s_ul_hash_t umis_per_gene;
		for (auto const &gene : this->_genes)
		{
			size_t umis_num = 0;
			for (auto const &umi : gene.second)
			{
				if (!umi.second.mark.match(this->_query_marks))
					continue;

				if (return_reads)
				{
					umis_num += umi.second.read_count;
				}
				else
				{
					umis_num++;
				}
			}

			if (umis_num == 0)
				continue;

			umis_per_gene.emplace(gene.first, umis_num);
		}

		return umis_per_gene;
	}

	Cell::ss_ul_hash_t Cell::requested_reads_per_umi_per_gene() const
	{
		ss_ul_hash_t reads_per_umi_per_gene;
		for (auto const &gene : this->_genes)
		{
			s_ul_hash_t reads_per_umi;
			for (auto const &umi : gene.second)
			{
				if (!umi.second.mark.match(this->_query_marks))
					continue;

				reads_per_umi.emplace(umi.first, umi.second.read_count);
			}

			if (reads_per_umi.empty())
				continue;

			reads_per_umi_per_gene.emplace(gene.first, reads_per_umi);
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

	const std::string &Cell::barcode() const
	{
		return this->_barcode;
	}

	size_t Cell::umis_number() const
	{
		return this->_umis_number;
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

	void Cell::update_requested_size()
	{
		this->_requested_genes_num = 0;
		this->_requested_umis_num = 0;
		for (auto const &gene : this->_genes)
		{
			bool has_umis = false;
			for (auto const &umi : gene.second)
			{
				if (umi.second.mark.match(this->_query_marks))
				{
					has_umis = true;
					this->_requested_umis_num++;
				}
			}

			if (has_umis)
			{
				this->_requested_genes_num++;
			}
		}
	}

	const Cell::umi_map_t& Cell::at(const std::string &gene) const
	{
		return this->_genes.at(gene);
	}
}