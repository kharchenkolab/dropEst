#include "RefGenesContainer.h"

#include <Tools/GeneInfo.h>
#include <Tools/Logs.h>

#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <string.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


namespace Tools
{
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

			auto iter = this->_genes_intervals.insert(std::make_pair(info.chr_name(), IntervalsContainer<GeneInfo>(true, 1)));
			iter.first->second.add_interval(info.start_pos(), info.end_pos(), info);
		}

		for (auto &intervals : _genes_intervals)
		{
			intervals.second.set_initialized();
		}
	}

	RefGenesContainer::gene_names_set_t RefGenesContainer::accumulate_genes(const genes_set_t &genes) const //Some genes are annotated only in pairs.
	{
		gene_names_set_t res;
		if (genes.empty())
			return res;

		const std::string &chr_name = genes.begin()->chr_name();
		for (auto const &gene : genes)
		{
			res.insert(gene.name());

			if (chr_name != gene.chr_name())
				throw std::runtime_error("Can't use genes from different chromosomes: " + genes.begin()->chr_name() +
				                         ", " + gene.chr_name());
		}

		return res;
	}

	RefGenesContainer::gene_names_set_t RefGenesContainer::get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const
	{
		if (end_pos < start_pos)
			return gene_names_set_t();

		auto current_intervals_it = this->_genes_intervals.find(chr_name);
		if (current_intervals_it == this->_genes_intervals.end())
			throw ChrNotFoundException(chr_name);

		auto const &current_intervals = current_intervals_it->second;
		auto intercepted_genes = current_intervals.get_intervals(start_pos, end_pos);

		genes_set_t res_genes;
		pos_t read_len = end_pos - start_pos;
		for (auto const &gene : intercepted_genes)
		{
			double intercept_len = std::min(gene.end_pos(), end_pos) - std::max(gene.start_pos(), start_pos);
			if (intercept_len / std::min(gene.length(), read_len) > RefGenesContainer::read_intersection_significant_part)
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

		if (id.empty())
		{
			if (name.empty())
				throw std::runtime_error("GTF record doesn't contain either gene name or id:\n" + record);

			id = name;
		}

		size_t start_pos = strtoul(columns[3].c_str(), NULL, 10) - 1;
		size_t end_pos = strtoul(columns[4].c_str(), NULL, 10);

		return GeneInfo(columns[0], id, name, start_pos, end_pos); // TODO: add transcript
	}

	GeneInfo RefGenesContainer::parse_bed_record(const std::string &record)
	{
		GeneInfo result;
		auto first_char_index = record.find_first_not_of("\t ");
		if (first_char_index == std::string::npos || record[first_char_index] == '#')
			return result;

		std::vector<std::string> columns(RefGenesContainer::split(record));
		if (columns.size() < 4)
			throw std::runtime_error("Bed record is too short:\n" + record);

		size_t start_pos = strtoul(columns[1].c_str(), NULL, 10);
		size_t end_pos = strtoul(columns[2].c_str(), NULL, 10);

		return GeneInfo(columns[0], columns[3], "", start_pos, end_pos);
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

	bool RefGenesContainer::is_empty() const
	{
		return this->_is_empty;
	}
}
