#include "RefGenesContainer.h"

#include <Tools/GeneInfo.h>
#include <Tools/Logs.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <string.h>


namespace Tools
{
	RefGenesContainer::RefGenesContainer(const std::string &gtf_filename)
	{
		this->init_from_gtf(gtf_filename);
	}

	void RefGenesContainer::init_from_gtf(const std::string &gtf_filename)
	{
		std::vector<std::unordered_map<std::string, genes_list_t>> genes_groups;

		std::string record;
		std::ifstream gtf_in(gtf_filename);
		while (std::getline(gtf_in, record))
		{
			GeneInfo info;
			try
			{
				info = RefGenesContainer::parse_gtf_record(record);
			}
			catch (std::runtime_error err)
			{
				L_ERR << err.what();
				continue;
			}

			if (!info.is_valid())
				continue;

			if (genes_groups.size() <= info.chr_num())
			{
				genes_groups.resize(info.chr_num() + 1);
			}
			RefGenesContainer::add_gene(info, genes_groups[info.chr_num()][info.id()]);
		}


		this->_genes_intervals.resize(genes_groups.size());

		for (size_t chr_id = 0; chr_id < genes_groups.size(); ++chr_id)
		{
			auto &cur_groups = genes_groups[chr_id];
			if (cur_groups.empty())
				continue;

			genes_vec_t genes;
			auto back_inserter = std::back_inserter(genes);
			for (auto &gene_group : cur_groups)
			{
				std::move(gene_group.second.begin(), gene_group.second.end(), back_inserter);
			}

			this->_genes_intervals[chr_id] = this->filter_genes(genes);
		}
	}

	void RefGenesContainer::add_gene(GeneInfo &gene, genes_list_t &genes)
	{
		auto cur_iterator = genes.begin();

		while (cur_iterator != genes.end() && !gene.is_intercept(*cur_iterator))
		{
			if (cur_iterator->start_pos() > gene.end_pos())
			{
				genes.insert(cur_iterator, gene);
				return;
			}
			++cur_iterator;
		}

		if (cur_iterator == genes.end())
		{
			genes.push_back(gene);
			return;
		}

		auto end_iterator = cur_iterator;
		end_iterator++;

		while (end_iterator != genes.end() && gene.is_intercept(*end_iterator))
		{
			gene.merge(*end_iterator);
			++end_iterator;
		}

		cur_iterator->merge(gene);
		cur_iterator++;
		genes.erase(cur_iterator, end_iterator);
	}

	RefGenesContainer::intervals_vec_t RefGenesContainer::filter_genes(const genes_vec_t &genes)
	{
		const int min_interval_len = 5;

		if (genes.empty())
			return intervals_vec_t();

		auto events = RefGenesContainer::genes_to_events(genes);

		intervals_vec_t result;

		pos_t start_pos = 0, end_pos;
		genes_set cur_genes;
		for (auto const &event : events)
		{
			end_pos = event.first;
			if (!cur_genes.empty() && end_pos - start_pos >= min_interval_len)
			{
				result.push_back(Interval{start_pos, end_pos, cur_genes});
				this->save_gene_names(cur_genes);
			}

			if (event.second->start_pos() == event.first)
			{
				cur_genes.insert(*event.second);
			}
			else
			{
				cur_genes.erase(*event.second);
			}

			start_pos = end_pos;
		}

		return result;
	}

	GeneInfo RefGenesContainer::accumulate_genes(const genes_set &genes) const
	{
		if (genes.empty())
			return GeneInfo();

		auto gene_it = genes.begin();
		std::string id = gene_it->id();
		if (genes.size() > 1 && this->single_gene_names.find(id) != this->single_gene_names.end())
			return GeneInfo();

		GeneInfo::num_t chr_num = gene_it->chr_num();
		pos_t end_pos = gene_it->end_pos();
		while (++gene_it != genes.end())
		{
			if (this->single_gene_names.find(gene_it->id()) != this->single_gene_names.end())
				return GeneInfo();

			id += "," + gene_it->id();

			if (chr_num != gene_it->chr_num())
				throw std::runtime_error("Can't use genes from different chromosomes: " + genes.begin()->chr_name() +
						                         ", " + gene_it->chr_name());

			end_pos = gene_it->end_pos();
		}

		return GeneInfo(genes.begin()->chr_name(), id, genes.begin()->start_pos(), end_pos);
	}

	GeneInfo RefGenesContainer::get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const
	{
		const double significant_part = 0.5;
		if (end_pos < start_pos)
			return GeneInfo();

		GeneInfo::num_t chr_num = GeneInfo::parse_chr_name(chr_name);
		const intervals_vec_t &current_intervals = this->_genes_intervals[chr_num];
		auto intercept_it = std::lower_bound(current_intervals.begin(), current_intervals.end(), start_pos,
											 [](const Interval &interval, pos_t pos){ return interval.end_pos <= pos;});

		if (intercept_it == current_intervals.end() || intercept_it->start_pos > end_pos)
			return GeneInfo();

		genes_set intercepted_genes;
		while (intercept_it != current_intervals.end() && intercept_it->start_pos <= end_pos)
		{
			intercepted_genes.insert(intercept_it->genes.begin(), intercept_it->genes.end());
			++intercept_it;
		}

		genes_set res_genes;
		pos_t read_len = end_pos - start_pos;
		for (auto const &gene : intercepted_genes)
		{
			double intercept_len = std::min(gene.end_pos(), end_pos) - std::max(gene.start_pos(), start_pos);
			if (intercept_len / gene.size() >= significant_part || intercept_len / read_len >= significant_part)
			{
				res_genes.insert(gene);
			}
		}

		return RefGenesContainer::accumulate_genes(res_genes);
	}

	GeneInfo RefGenesContainer::parse_gtf_record(const std::string &record)
	{
		GeneInfo result;
		std::istringstream istr(record);

		std::string column;
		std::vector<std::string> columns;
		while (istr >> column)
		{
			columns.push_back(column);
		}

		if (columns.size() < 9)
			throw std::runtime_error("Can't parse record: \n" + record);

		if (columns[0] == "." || columns[3] == "." || columns[4] == "." || columns.size() == 9)
			return result;

		std::string id = "";
		for (size_t attrib_ind = 8; attrib_ind < columns.size(); ++attrib_ind)
		{
			std::string key = columns[attrib_ind];
			std::string value = columns[attrib_ind + 1];

			if (key == "gene_id")
			{
				id = value.substr(1, value.length() - 3);
				break;
			}
		}

		return GeneInfo(columns[0], id, strtoul(columns[3].c_str(), NULL, 10), strtoul(columns[4].c_str(), NULL, 10));
	}

	RefGenesContainer::gene_event_t RefGenesContainer::genes_to_events(const genes_vec_t &genes)
	{
		size_t chr_num = 0;
		if (!genes.empty())
		{
			chr_num = genes.front().chr_num();
		}

		gene_event_t events;
		for (auto const &gene : genes)
		{
			if (chr_num != gene.chr_num())
				throw std::runtime_error("Genes from different chromosomes in one group: " + gene.chr_name() + ", " +
				                         gene.id() + " and " + genes.front().chr_name() + ", " + genes.front().id());

			events.emplace(gene.start_pos(), &gene);
			events.emplace(gene.end_pos(), &gene);
		}

		return events;
	}

	void RefGenesContainer::save_gene_names(const RefGenesContainer::genes_set &genes_in_interval)
	{
		if (genes_in_interval.size() != 1)
			return;

		this->single_gene_names.insert(genes_in_interval.begin()->id());
	}
}
