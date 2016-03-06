#pragma once

#include <list>
#include <set>
#include <string>
#include <vector>

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

		class GeneInfo
		{
			friend struct TestTools::testGtf;
			friend struct TestTools::testGeneMerge;
		public:
			typedef unsigned int num_t;
		private:
			std::string _chr_name = "";
			std::string _id = "";
			pos_t _start_pos;
			pos_t _end_pos;
			num_t _chr_num = 0;

		public:
			GeneInfo() = default;
			GeneInfo(const std::string &chr_name, std::string id, pos_t start_pos, pos_t end_pos);

			std::string chr_name() const;
			num_t chr_num() const;
			std::string id() const;
			pos_t start_pos() const;
			pos_t end_pos() const;

			bool is_valid() const;
			void merge(const GeneInfo &other);
			bool is_intercept(const GeneInfo &other) const;

			bool operator<(const GeneInfo &other) const;

			static num_t parse_chr_name(const std::string &chr_name);
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
		typedef std::vector<GeneInfo> genes_vec_t;
		typedef std::vector<Interval> intervals_vec_t;

	private:
		std::vector<intervals_vec_t> _genes_intervals;


	private:
		void init_from_gtf(const std::string &gtf_filename);
		static void add_gene(GeneInfo &gene, genes_list_t &genes);
		static GeneInfo parse_gtf_record(const std::string &record);
		static GeneInfo accumulate_genes(genes_set genes);

	public:
		RefGenesContainer(const std::string &gtf_filename);
		GeneInfo get_gene_info(const std::string &chr_name, pos_t start_pos, pos_t end_pos) const;

		intervals_vec_t filter_genes(const genes_vec_t &genes);
	};
}