#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Gene.h"
#include "Stats.h"
#include "StringIndexer.h"

namespace Estimation
{
	class Cell
	{
	public:
		typedef std::map<StringIndexer::index_t, Gene> genes_t;
		typedef std::unordered_map<std::string, std::string> s_s_hash_t;
		typedef std::unordered_map<std::string, size_t> s_ul_hash_t;
		typedef std::unordered_map<std::string, s_ul_hash_t> ss_ul_hash_t;

	private:
		std::unique_ptr<char[]> _barcode;
		const size_t _min_genes_to_be_real;

		bool _is_merged;
		bool _is_excluded;
		size_t _requested_genes_num;
		size_t _requested_umis_num;

		genes_t _genes;
		Stats _stats;

		StringIndexer *_gene_indexer;
		StringIndexer *_umi_indexer;

	public:
		bool is_merged() const;
		bool is_excluded() const;
		bool is_real() const;
		std::string barcode() const;
		const char* barcode_c() const;
		size_t umis_number() const;
		size_t requested_genes_num() const;
		size_t requested_umis_num() const;

		const Stats &stats() const;
		Stats &stats();
		const genes_t& genes() const;
		s_ul_hash_t requested_umis_per_gene(const UMI::Mark::query_t &query_marks, bool return_reads) const;
		ss_ul_hash_t requested_reads_per_umi_per_gene(const UMI::Mark::query_t &query_marks) const;

		size_t size() const;
		const Gene& at(const std::string &gene) const;

		void add_umi(const std::string &gene, const std::string &umi, const UMI::Mark &umi_mark);
		void set_merged();
		void set_excluded();
		void merge(const Cell &source);
		void merge_umis(StringIndexer::index_t gene, const s_s_hash_t &merge_targets);
		void update_requested_size(const UMI::Mark::query_t &query_marks);

		Cell(const std::string &barcode, size_t min_genes_to_be_real, StringIndexer *gene_indexer,
		     StringIndexer *umi_indexer);
	};
}
