#pragma once

#include "GtfRecord.h"
#include <Tools/IntervalsContainer.h>

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace TestTools
{
	struct testGtf;
	struct testInitGtf;
	struct testGeneMerge;
	struct testParseBed;
}

namespace Tools
{
	class RefGenesContainer
	{
		friend struct TestTools::testGtf;
		friend struct TestTools::testInitGtf;
		friend struct TestTools::testGeneMerge;
		friend struct TestTools::testParseBed;

	public:
		typedef unsigned long pos_t;

		class ChrNotFoundException : public std::runtime_error
		{
		public:
			const std::string chr_name;
			ChrNotFoundException(const std::string &chr_name)
				: std::runtime_error("Can't find chromosome " + chr_name)
				, chr_name(chr_name)
			{}
		};

		class QueryResult
		{
		public:
			std::string gene_name;
			GtfRecord::RecordType type;

			bool operator<(const QueryResult &other) const;
			QueryResult(const std::string &gene_name = "", GtfRecord::RecordType type = GtfRecord::NONE);
		};

		typedef std::set<QueryResult> query_results_t;

	private:
		typedef std::set<GtfRecord> genes_set_t;
		typedef std::unordered_map<std::string, IntervalsContainer<std::string>> transcript_intervals_map_t;
		typedef std::unordered_map<std::string, std::unordered_map<std::string, IntervalsContainer<GtfRecord::RecordType>>> exon_intervals_map_t;
		typedef std::unordered_map<std::string, std::unordered_map<std::string, Interval>> transcript_positions_map_t;

	private:
		bool _is_empty;
		std::string _file_format;

		bool _use_introns_from_gtf;
		bool _gtf_has_transcripts;

		transcript_intervals_map_t _transcript_intervals; // chr -> IntervalsContainer<transcript>
		exon_intervals_map_t _exons_by_transcripts; //chr -> transcript -> IntervalsContainer<PositionType>
		transcript_positions_map_t _transcript_positions; //chr -> transcript -> [start, end)
		std::unordered_map<std::string, std::string> _genes_by_transcripts;


	private:
		void init(const std::string &genes_filename);
		void save_transcript(GtfRecord record);
		GtfRecord parse_gtf_record(const std::string &record);

		static GtfRecord parse_bed_record(const std::string &record);
		static std::vector<std::string> split(const std::string &record);

	public:
		RefGenesContainer();
		RefGenesContainer(const std::string &genes_filename);

		/// Get genes, intersected requested interval
		/// \param chr_name chromosome name
		/// \param start_pos start position in 0-based coordinate system (inclusive)
		/// \param end_pos end position in 0-based coordinate system (exclusive)
		/// \return
		query_results_t get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const;
		bool is_empty() const;
		bool has_introns() const;
	};
}
