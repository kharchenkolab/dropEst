#include "RefGenesContainer.h"

#include "GtfRecord.h"
#include <Tools/Logs.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string.h>
#include <unordered_map>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace Tools
{
	namespace GeneAnnotation
	{
		RefGenesContainer::RefGenesContainer()
				: _is_empty(true)
				, _use_introns_from_gtf(false)
				, _gtf_has_transcripts(true)
		{}

		RefGenesContainer::RefGenesContainer(const std::string &genes_filename)
				: _is_empty(false)
				, _use_introns_from_gtf(false)
				, _gtf_has_transcripts(true)
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
			std::ifstream gtf_in(genes_filename);
			if (gtf_in.fail())
				throw std::runtime_error("Can't open GTF file: '" + genes_filename + "'");

			boost::iostreams::filtering_istream gz_fs;
			if (genes_filename.substr(genes_filename.length() - 3) == ".gz")
			{
				gz_fs.push(boost::iostreams::gzip_decompressor());
			}
			gz_fs.push(gtf_in);

			std::string line;
			while (std::getline(gz_fs, line))
			{
				GtfRecord record;
				try
				{
					if (this->_file_format == "gtf")
					{
						record = RefGenesContainer::parse_gtf_record(line);
					} else
					{
						record = RefGenesContainer::parse_bed_record(line);
					}
				}
				catch (std::runtime_error err)
				{
					L_ERR << err.what();
					continue;
				}

				if (!record.is_valid())
					continue;

				this->save_transcript(record);
			}

			for (auto const &chr : this->_transcript_positions)
			{
				auto &intervals = this->_transcript_intervals.emplace(chr.first, true).first->second;
				for (auto const &transcript : chr.second)
				{
					intervals.add_interval(transcript.second.start_pos(), transcript.second.end_pos(),
					                       transcript.first);
					this->_exons_by_transcripts.at(chr.first).at(transcript.first).set_initialized();
				}
				intervals.set_initialized();
			}
		}

		void RefGenesContainer::save_transcript(GtfRecord record)
		{
			auto transcript_iter = this->_transcript_positions[record.chr_name()].insert(
					std::make_pair(record.transcript_id(), record));
			transcript_iter.first->second.merge(record);

			auto exon_iter = this->_exons_by_transcripts[record.chr_name()].emplace(record.transcript_id(),
			                                                                        IntervalsContainer<GtfRecord::RecordType>(
					                                                                        false));
			exon_iter.first->second.add_interval(record.start_pos(), record.end_pos(), record.type());

			auto gene_iter = this->_genes_by_transcripts.emplace(record.transcript_id(), record.gene_name());
			if (!gene_iter.second && gene_iter.first->second != record.gene_name())
				throw std::runtime_error(
						"Different gene names (" + record.gene_name() + ", " + gene_iter.first->second +
						") for the same transcript (" + record.transcript_id() + ")");
		}

		GtfRecord RefGenesContainer::parse_gtf_record(const std::string &record)
		{
			GtfRecord result;
			if (record.at(0) == '#')
				return result;

			std::vector<std::string> columns(RefGenesContainer::split(record));

			if (columns.size() < 9)
				throw std::runtime_error("Can't parse record: \n" + record);

			if (columns[0] == "." || columns[3] == "." || columns[4] == "." || columns.size() == 9)
				return result;

			GtfRecord::RecordType type;
			if (columns[2] == "exon")
			{
				type = GtfRecord::EXON;
			} else if (columns[2] == "intron")
			{
				type = GtfRecord::INTRON;
				this->_use_introns_from_gtf = true;
			} else
			{
				return result;
			}

			std::string id, name, transcript;
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
				if (key == "transcript_id")
				{
					transcript = value.substr(1, value.length() - 3);
				}
			}

			if (transcript.empty())
			{
				this->_gtf_has_transcripts = false;
			}

			if (id.empty())
			{
				if (name.empty())
					throw std::runtime_error("GTF record doesn't contain either gene name or id:\n" + record);

				id = name;
			}

			size_t start_pos = strtoul(columns[3].c_str(), NULL, 10) - 1;
			size_t end_pos = strtoul(columns[4].c_str(), NULL, 10);

			return GtfRecord(columns[0], id, name, start_pos, end_pos, type, transcript);
		}

		RefGenesContainer::query_results_t
		RefGenesContainer::get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const
		{
			if (end_pos < start_pos)
				return query_results_t();

			auto current_intervals_it = this->_transcript_intervals.find(chr_name);
			if (current_intervals_it == this->_transcript_intervals.end())
				throw ChrNotFoundException(chr_name);

			auto const &current_intervals = current_intervals_it->second;

			query_results_t results;
			auto intersected_transcripts = current_intervals.get_intervals(start_pos, end_pos);
			for (const std::string &transcript : intersected_transcripts)
			{
				auto position_types = this->_exons_by_transcripts.at(chr_name).at(transcript).get_intervals(start_pos,
				                                                                                            end_pos); // TODO: it works only for single-nucleotide queries. After CIGAR parsing it should be rewritten.
				const std::string &gene_name = this->_genes_by_transcripts.at(transcript);
				if (position_types.empty() && !this->_use_introns_from_gtf)
				{
					results.emplace(gene_name, GtfRecord::INTRON);
					continue;
				}
				for (auto &type : position_types)
				{
					results.emplace(gene_name, type);
				}
			}

			return results;
		}

		GtfRecord RefGenesContainer::parse_bed_record(const std::string &record)
		{
			GtfRecord result;
			auto first_char_index = record.find_first_not_of("\t ");
			if (first_char_index == std::string::npos || record[first_char_index] == '#')
				return result;

			std::vector<std::string> columns(RefGenesContainer::split(record));
			if (columns.size() < 4)
				throw std::runtime_error("Bed record is too short:\n" + record);

			size_t start_pos = strtoul(columns[1].c_str(), NULL, 10);
			size_t end_pos = strtoul(columns[2].c_str(), NULL, 10);

			return GtfRecord(columns[0], columns[3], "", start_pos, end_pos, GtfRecord::EXON);
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

		bool RefGenesContainer::has_introns() const
		{
			return this->_gtf_has_transcripts || this->_use_introns_from_gtf;
		}

		bool RefGenesContainer::QueryResult::operator<(const QueryResult &other) const
		{
			if (this->type == other.type)
				return this->gene_name < other.gene_name;

			return this->type < other.type;
		}

		RefGenesContainer::QueryResult::QueryResult(const std::string &gene_name, GtfRecord::RecordType type)
				: gene_name(gene_name)
				, type(type)
		{}
	}
}
