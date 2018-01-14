#pragma once

#include "GtfRecord.h"
#include "IntervalsContainer.h"

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
	namespace GeneAnnotation
	{
		class RefGenesContainer
		{
			friend struct TestTools::testGtf;
			friend struct TestTools::testInitGtf;
			friend struct TestTools::testGeneMerge;
			friend struct TestTools::testParseBed;

		public:
			using pos_t = unsigned long;

			class ChrNotFoundException : public std::runtime_error
			{
			public:
				const std::string chr_name;

				explicit ChrNotFoundException(const std::string &chr_name)
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

				explicit QueryResult(const std::string &gene_name = "", GtfRecord::RecordType type = GtfRecord::NONE);
			};

			using query_results_t = std::set<QueryResult>;

		private:
			using transcript_intervals_map_t = std::unordered_map<std::string, IntervalsContainer<std::string>>;
			using exon_intervals_map_t = std::unordered_map<std::string, std::unordered_map<std::string, IntervalsContainer<GtfRecord::RecordType>>>;
			using transcript_positions_map_t = std::unordered_map<std::string, std::unordered_map<std::string, Interval>>;

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

			explicit RefGenesContainer(const std::string &genes_filename);

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
}
