#pragma once

#include "GeneInfo.h"

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
}

namespace Tools
{
	class RefGenesContainer
	{
		friend struct TestTools::testGtf;
		friend struct TestTools::testInitGtf;
		friend struct TestTools::testGeneMerge;

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

	private:
		typedef std::set<GeneInfo> genes_set;
		struct Interval
		{
			pos_t start_pos;
			pos_t end_pos;
			genes_set genes;
		};

		typedef std::list<GeneInfo> genes_list_t;
		typedef std::multimap<pos_t, const GeneInfo*> gene_event_t;
		typedef std::vector<GeneInfo> genes_vec_t;
		typedef std::vector<Interval> intervals_vec_t;
		typedef std::unordered_map<std::string, intervals_vec_t> intervals_map_t;

	private:
		static const int min_interval_len;
		static const double read_intersection_significant_part;

		intervals_map_t _genes_intervals;
		std::unordered_set<std::string> _single_gene_names;
		bool _is_empty;


	private:
		void init_from_gtf(const std::string &gtf_filename);
		static void add_gene(GeneInfo &gene, genes_list_t &genes);
		static GeneInfo parse_gtf_record(const std::string &record);
		GeneInfo accumulate_genes(const genes_set &genes) const;

		static gene_event_t genes_to_events(const genes_vec_t &genes);

		intervals_vec_t filter_genes(const genes_vec_t &genes);

	public:
		RefGenesContainer();
		RefGenesContainer(const std::string &gtf_filename);
		GeneInfo get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const;
		bool is_empty() const;

		void save_gene_names(const genes_set &genes_in_interval);
	};
}
