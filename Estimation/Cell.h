#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "Gene.h"

namespace Estimation
{
	class Cell
	{
	public:
		typedef std::map<std::string, Gene> genes_t;
		typedef std::unordered_map<std::string, std::string> s_s_hash_t;
		typedef std::unordered_map<std::string, size_t> s_ul_hash_t;
		typedef std::unordered_map<std::string, s_ul_hash_t> ss_ul_hash_t;

	private:
		const std::string _barcode;
		const size_t _min_genes_to_be_real;
		const UMI::Mark::query_t _query_marks;


		bool _is_merged;
		bool _is_excluded;
		size_t _umis_number;
		size_t _requested_genes_num;
		size_t _requested_umis_num;

		genes_t _genes;

	public:
		bool is_merged() const;
		bool is_excluded() const;
		bool is_real() const;
		const std::string& barcode() const;
		size_t umis_number() const;
		size_t requested_genes_num() const;
		size_t requested_umis_num() const;

		void add_umi(const std::string &gene, const std::string &umi, const UMI::Mark &umi_mark);
		void set_merged();
		void set_excluded();
		void merge(const Cell &source);
		void merge_umis(const std::string &gene, const s_s_hash_t &merge_targets);
		void update_requested_size();

		const genes_t& genes() const;
		s_ul_hash_t requested_umis_per_gene(bool return_reads) const;
		ss_ul_hash_t requested_reads_per_umi_per_gene() const;

		size_t size() const;
		const Gene& at(const std::string &gene) const;

		Cell(const std::string &barcode, size_t min_genes_to_be_real, const std::vector<UMI::Mark> &query_marks);
	};
}
