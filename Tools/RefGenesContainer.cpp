#include "RefGenesContainer.h"

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
		time_t t_parse = 0, t_add = 0;
		while (std::getline(gtf_in, record))
		{
			GeneInfo info;
			try
			{
				time_t t0 = clock();
				info = RefGenesContainer::parse_gtf_record(record);
				t_parse += clock() - t0;
			}
			catch (std::runtime_error err)
			{
				L_ERR << err.what();
				continue;
			}

			if (!info.is_valid())
				continue;
			time_t t1 = clock();
			if (genes_groups.size() <= info.chr_num())
			{
				genes_groups.resize(info.chr_num() + 1);
			}
			RefGenesContainer::add_gene(info, genes_groups[info.chr_num()][info.id()]);
			t_add += clock() - t1;
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

			this->_genes_intervals[chr_id] = RefGenesContainer::filter_genes(genes);
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
		if (genes.empty())
			return intervals_vec_t();

		intervals_vec_t result;
		std::multimap<pos_t, const GeneInfo*> events;
		for (auto const &gene : genes)
		{
			events.emplace(gene.start_pos(), &gene);
			events.emplace(gene.end_pos(), &gene);
		}

		pos_t start_pos = 0, end_pos;
		genes_set cur_genes;
		for (auto const &event : events)
		{
			end_pos = event.first;
			if (!cur_genes.empty() && start_pos != end_pos)
			{
				result.push_back(Interval{start_pos, end_pos, cur_genes});
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

	RefGenesContainer::GeneInfo RefGenesContainer::accumulate_genes(RefGenesContainer::genes_set genes)
	{
		if (genes.empty())
			return GeneInfo();

		auto gene_it = genes.begin();
		std::string id = gene_it->id();
		GeneInfo::num_t chr_num = gene_it->chr_num();
		pos_t end_pos = gene_it->end_pos();
		while (++gene_it != genes.end())
		{
			id += "," + gene_it->id();

			if (chr_num != gene_it->chr_num())
			{
				throw std::runtime_error("Can't use genes from different chromosomes: " + genes.begin()->chr_name() + ", " + gene_it->chr_name());
			}
			end_pos = gene_it->end_pos();
		}

		return GeneInfo(genes.begin()->chr_name(), id, genes.begin()->start_pos(), end_pos);
	}

	RefGenesContainer::GeneInfo RefGenesContainer::get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const
	{
		GeneInfo::num_t chr_num = GeneInfo::parse_chr_name(chr_name);
		const intervals_vec_t &current_intervals = this->_genes_intervals[chr_num];
		auto intercept_it = std::lower_bound(current_intervals.begin(), current_intervals.end(), start_pos,
											 [](const Interval &interval, pos_t pos){ return interval.end_pos <= pos;});

		if (intercept_it == current_intervals.end() || intercept_it->start_pos > end_pos)
			return GeneInfo();

		genes_set res_genes;
		while (intercept_it != current_intervals.end() && intercept_it->start_pos <= end_pos)
		{
			res_genes.insert(intercept_it->genes.begin(), intercept_it->genes.end());
			++intercept_it;
		}

		return RefGenesContainer::accumulate_genes(res_genes);
	}

	RefGenesContainer::GeneInfo RefGenesContainer::parse_gtf_record(const std::string &record)
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

	bool RefGenesContainer::GeneInfo::is_valid() const
	{
		return this->id() != "";
	}

	bool RefGenesContainer::GeneInfo::operator<(const GeneInfo &other) const
	{
		return this->id() < other.id();
	}

	bool RefGenesContainer::GeneInfo::is_intercept(const GeneInfo &other) const
	{
		return this->start_pos() <= other.end_pos() && this->end_pos() >= other.start_pos();
	}

	void RefGenesContainer::GeneInfo::merge(const GeneInfo &other)
	{
		this->_start_pos = std::min(this->_start_pos, other._start_pos);
		this->_end_pos= std::max(this->_end_pos, other._end_pos);
	}

	RefGenesContainer::GeneInfo::GeneInfo(const std::string &chr_name, std::string id, pos_t start_pos, pos_t end_pos)
		: _chr_name(chr_name)
		, _id(id)
		, _start_pos(start_pos)
		, _end_pos(end_pos)
	{
		this->_chr_num = GeneInfo::parse_chr_name(chr_name);
	}

	std::string RefGenesContainer::GeneInfo::chr_name() const
	{
		return this->_chr_name;
	}

	std::string RefGenesContainer::GeneInfo::id() const
	{
		return this->_id;
	}

	RefGenesContainer::pos_t RefGenesContainer::GeneInfo::start_pos() const
	{
		return this->_start_pos;
	}

	RefGenesContainer::pos_t RefGenesContainer::GeneInfo::end_pos() const
	{
		return this->_end_pos;
	}

	RefGenesContainer::GeneInfo::num_t RefGenesContainer::GeneInfo::chr_num() const
	{
		return this->_chr_num;
	}

	RefGenesContainer::GeneInfo::num_t RefGenesContainer::GeneInfo::parse_chr_name(const std::string &chr_name)
	{
		return chr_name.length() > 2 ? atoi(chr_name.substr(3).c_str()) : atoi(chr_name.c_str());
	}
}
