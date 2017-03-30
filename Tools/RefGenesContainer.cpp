#include "RefGenesContainer.h"

#include <Tools/GeneInfo.h>
#include <Tools/Logs.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <string.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


namespace Tools
{
	const int RefGenesContainer::min_interval_len = 5;
	const double RefGenesContainer::read_intersection_significant_part = 0.0; // TODO: it doesn't work because of single-nucleotide queries in ReadsParamsParser

	RefGenesContainer::RefGenesContainer()
		: _is_empty(true)
	{}

	RefGenesContainer::RefGenesContainer(const std::string &genes_filename)
		: _is_empty(false)
	{
		auto wrong_format_exception = std::runtime_error("Wrong genes file format: '" + genes_filename + "'");
		if (genes_filename.length() < 3)
			throw wrong_format_exception;

		this->_file_format = genes_filename.substr(genes_filename.length() - 3);
		if (this->_file_format == ".gz")
		{
			if (genes_filename.length() < 6)
				throw wrong_format_exception;

			this->_file_format = genes_filename.substr(genes_filename.length() - 6, 3);
		}

		if (this->_file_format != "bed" && this->_file_format != "gtf")
			throw wrong_format_exception;

		this->init(genes_filename);
	}

	void RefGenesContainer::init(const std::string &genes_filename)
	{
		std::unordered_map<std::string, std::unordered_map<std::string, genes_list_t>> genes_by_chr_by_id;

		std::string record;
		std::ifstream gtf_in(genes_filename);
		if (gtf_in.fail())
			throw std::runtime_error("Can't open GTF file: '" + genes_filename + "'");

		boost::iostreams::filtering_istream gz_fs;
		if (genes_filename.substr(genes_filename.length() - 3) == ".gz")
		{
			gz_fs.push(boost::iostreams::gzip_decompressor());
		}
		gz_fs.push(gtf_in);

		while (std::getline(gz_fs, record))
		{
			GeneInfo info;
			try
			{
				if (this->_file_format == "gtf")
				{
					info = RefGenesContainer::parse_gtf_record(record);
				}
				else
				{
					info = RefGenesContainer::parse_bed_record(record);
				}
			}
			catch (std::runtime_error err)
			{
				L_ERR << err.what();
				continue;
			}

			if (!info.is_valid())
				continue;

			RefGenesContainer::add_gene(info, genes_by_chr_by_id[info.chr_name()][info.id()]);
		}

		for (auto const &chr_genes : genes_by_chr_by_id)
		{
			if (chr_genes.second.empty())
				continue;

			genes_vec_t genes;
			auto back_inserter = std::back_inserter(genes);
			for (auto &gene_group : chr_genes.second)
			{
				std::move(gene_group.second.begin(), gene_group.second.end(), back_inserter);
			}

			this->_genes_intervals[chr_genes.first] = this->filter_genes(genes);
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

		auto events = RefGenesContainer::genes_to_events(genes);

		intervals_vec_t result; // Interval is a continuous part of a reference with the same genes composition

		pos_t start_pos = 0, end_pos;
		genes_set cur_genes;
		for (auto const &event : events)
		{
			end_pos = event.first;
			if (!cur_genes.empty() && end_pos - start_pos >= RefGenesContainer::min_interval_len)
			{
				result.push_back(Interval{start_pos, end_pos, cur_genes});
				this->save_gene_names(cur_genes);
			}

			if (event.second->start_pos() == event.first) // "Gene start" event
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

	GeneInfo RefGenesContainer::accumulate_genes(const genes_set &genes) const //Some genes are annotated only in pairs.
	{
		if (genes.empty())
			return GeneInfo();

		auto gene_it = genes.begin();
		std::string id = gene_it->id(), name = gene_it->name();

		if (genes.size() > 1 && this->_single_gene_names.find(id) != this->_single_gene_names.end())
			return GeneInfo();

		const std::string &chr_name = gene_it->chr_name();
		pos_t end_pos = gene_it->end_pos();
		while (++gene_it != genes.end())
		{
			if (this->_single_gene_names.find(gene_it->id()) != this->_single_gene_names.end())
				return GeneInfo();

			id += "," + gene_it->id();
			name += "," + gene_it->name();

			if (chr_name != gene_it->chr_name())
				throw std::runtime_error("Can't use genes from different chromosomes: " + genes.begin()->chr_name() +
						                         ", " + gene_it->chr_name());

			end_pos = gene_it->end_pos();
		}

		return GeneInfo(genes.begin()->chr_name(), id, name, genes.begin()->start_pos(), end_pos);
	}

	GeneInfo RefGenesContainer::get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const
	{
		if (end_pos < start_pos)
			return GeneInfo();

		auto current_intervals_it = this->_genes_intervals.find(chr_name);
		if (current_intervals_it == this->_genes_intervals.end())
			throw ChrNotFoundException(chr_name);

		auto const &current_intervals = current_intervals_it->second;
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
			if (intercept_len / std::min(gene.size(), read_len) > RefGenesContainer::read_intersection_significant_part)
			{
				res_genes.insert(gene);
			}
		}

		return RefGenesContainer::accumulate_genes(res_genes);
	}

	GeneInfo RefGenesContainer::parse_gtf_record(const std::string &record)
	{
		GeneInfo result;
		if (record.at(0) == '#')
			return result;

		std::vector<std::string> columns(RefGenesContainer::split(record));

		if (columns.size() < 9)
			throw std::runtime_error("Can't parse record: \n" + record);

		if (columns[0] == "." || columns[3] == "." || columns[4] == "." || columns.size() == 9)
			return result;

		std::string id, name;
		for (size_t attrib_ind = 8; attrib_ind < columns.size() - 1; ++attrib_ind)
		{
			std::string key = columns[attrib_ind];
			std::string value = columns[attrib_ind + 1];

			if (key == "gene_id")
			{
				id = value.substr(1, value.length() - 3);
			}

			if (key == "gene_name")
			{
				name = value.substr(1, value.length() - 3);
			}
		}

		return GeneInfo(columns[0], id, name, strtoul(columns[3].c_str(), NULL, 10) - 1, strtoul(columns[4].c_str(), NULL, 10));
	}

	GeneInfo RefGenesContainer::parse_bed_record(const std::string &record)
	{
		GeneInfo result;
		if (record.at(0) == '#')
			return result;

		std::vector<std::string> columns(RefGenesContainer::split(record));
		return GeneInfo(columns[0], columns[3], "", strtoul(columns[1].c_str(), NULL, 10), strtoul(columns[2].c_str(), NULL, 10));
	}

	std::vector<std::string> RefGenesContainer::split(const std::string &record)
	{
		std::istringstream istr(record);

		std::string column;
		std::vector<std::string> columns;
		while (istr >> column)
		{
			columns.push_back(column);
		}
		return columns;
	}

	RefGenesContainer::gene_events_t RefGenesContainer::genes_to_events(const genes_vec_t &genes)
	{
		std::string chr_name;
		if (!genes.empty())
		{
			chr_name = genes.front().chr_name();
		}

		gene_events_t events;
		for (auto const &gene : genes)
		{
			if (chr_name != gene.chr_name())
				throw std::runtime_error("Genes from different chromosomes in one group: " + gene.chr_name() + ", " +
				                         gene.id() + " and " + chr_name + ", " + genes.front().id());

			events.emplace(gene.start_pos(), &gene);
			events.emplace(gene.end_pos(), &gene);
		}

		return events;
	}

	void RefGenesContainer::save_gene_names(const RefGenesContainer::genes_set &genes_in_interval)
	{
		if (genes_in_interval.size() != 1)
			return;

		this->_single_gene_names.insert(genes_in_interval.begin()->id());
	}

	bool RefGenesContainer::is_empty() const
	{
		return this->_is_empty;
	}
}
