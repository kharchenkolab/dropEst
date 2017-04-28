#pragma once

#include "GeneInfo.h"
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
		typedef std::set<std::string> gene_names_set_t;

		class ChrNotFoundException : public std::runtime_error
		{
		public:
			const std::string chr_name;
			ChrNotFoundException(const std::string &chr_name)
				: std::runtime_error("Can't find chromosome " + chr_name)
				, chr_name(chr_name)
			{}
		};

	private:
		typedef std::set<GeneInfo> genes_set_t;
		typedef std::unordered_map<std::string, IntervalsContainer<GeneInfo>> intervals_map_t;

	private:
		static const double read_intersection_significant_part; // TODO: move to parameter

		intervals_map_t _genes_intervals;
		bool _is_empty;
		std::string _file_format;


	private:
		void init(const std::string &genes_filename);
		static GeneInfo parse_gtf_record(const std::string &record);
		static GeneInfo parse_bed_record(const std::string &record);
		static std::vector<std::string> split(const std::string &record);

		gene_names_set_t accumulate_genes(const genes_set_t &genes) const;

	public:
		RefGenesContainer();
		RefGenesContainer(const std::string &genes_filename);

		/// Get genes, intersected requested interval
		/// \param chr_name chromosome name
		/// \param start_pos start position in 0-based coordinate system (inclusive)
		/// \param end_pos end position in 0-based coordinate system (exclusive)
		/// \return
		gene_names_set_t get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const;
		bool is_empty() const;
	};
}
