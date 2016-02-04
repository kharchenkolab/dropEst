#pragma once

#include "Stats.h"

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <vector>

#include <bam.h>

class GenesContainer
{
private:
	typedef boost::unordered_set<size_t> i_set_t;
	typedef boost::unordered_map<int, int> i_i_hash_t;
	typedef boost::unordered_map<std::string, i_i_hash_t> s_ii_hash_t;
	typedef boost::unordered_map<size_t, i_set_t> ISIHM;

public:
	struct IndexedCount
	{
		IndexedCount(size_t index, size_t count)
				: index(index)
				, count(count)
		{}

		size_t index;
		size_t count;

		static bool counts_comp(const IndexedCount &ic1, const IndexedCount &ic2)
		{
			return ic1.count < ic2.count;
		}
	};

	typedef boost::unordered_map<std::string, int> s_i_hash_t;
	typedef boost::unordered_map<std::string, s_i_hash_t> genes_t;

	typedef std::vector<IndexedCount> i_counter_t;
	typedef std::vector<int> ints_t;
	typedef std::vector<size_t> ids_t;
	typedef std::vector<std::string> names_t;

private:
	size_t top_print_size;

	std::vector<genes_t> _cells_genes; //cell_id -> gen_name -> umi -> count
	names_t _genes_names;
	ids_t _filtered_cells;

	i_counter_t _cells_genes_counts_sorted;

	Stats _stats;

private:
	void parse_bam_file(const std::string &bam_file_name, size_t read_prefix_length, s_i_hash_t &cells_ids, s_ii_hash_t &umig_cbs);
	bool parse_read_name(const bam1_t *align_info, std::string &read_name, std::string &cell_barcode, std::string &umi) const;

	void merge_genes(const s_ii_hash_t &umig_cbs, double min_merge_fraction,
					  int min_genes_after_merge, int max_merge_edit_distance);
	size_t get_umig_top(const ints_t &cb_reassigned, const IndexedCount &gene_count,
						const s_ii_hash_t &umigs_cells_counts, i_i_hash_t &umig_top) const;

	i_counter_t count_cells_genes(int min_genes_before_merge, bool logs = true) const;

	std::string get_iseq_verbose(bam1_t *align_info, size_t read_prefix_length) const;
	std::string get_cb_count_top_verbose(const i_counter_t &cells_genes_counts) const;
public:

	GenesContainer(const std::vector<std::string> &files, bool merge_tags, size_t read_prefix_length,
				   double min_merge_fraction, int min_genes_before_merge, int min_genes_after_merge,
				   int max_merge_edit_distance, size_t top_print_size);

	const ids_t& filtered_cells() const;
	const Stats& stats() const;
	const i_counter_t& cells_genes_counts_sorted() const;
	const genes_t& cell_genes(size_t index) const;
	const std::string& gene_name(size_t index) const;
};