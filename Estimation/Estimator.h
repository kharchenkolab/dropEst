#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <bam.h>
#include <boost/functional/hash.hpp>
#include <boost/property_tree/ptree.hpp>

#include "IndropResults.h"

class Estimator
{
private:
	typedef std::unordered_map<std::string, int, boost::hash<std::string> > SIHM;
	typedef std::unordered_map<std::string, SIHM, boost::hash<std::string> > SHHM;

	typedef std::unordered_map<int, int> IIHM;
	typedef std::unordered_map<std::string, IIHM> SIIHM;
	typedef std::unordered_map<int, std::unordered_set<int> > ISIHM;

	typedef std::vector<std::pair<int, int>> i_counter_t;
	typedef std::vector<std::pair<std::string, int>> s_counter_t;
	typedef std::vector<int> ints_t;
	typedef std::vector<std::string> names_t;

private:
	const size_t top_print_size = 10;

	size_t read_prefix_length;
	double min_merge_fraction;
	int max_merge_edit_distance;

	int min_genes;
	int low_genes;

	SIHM cb_ids, nonexone_chrs, exone_chrs;

	std::vector<SHHM> cb_genes;
	std::vector<std::string> cb_names;
	SIIHM umig_cbs;

public:
	Estimator(const std::vector<std::string> &files, const boost::property_tree::ptree &config);

	void run(bool merge_tags, bool text_output, const std::string &output_name);

private:
	void parse_bam_file(const std::string &bam_file_name);
	void merge_genes_with_tags(const i_counter_t &cb_genes_count_by_id, ints_t &merge_n, ints_t &unmerged_cbs); //TODO!!!

	void merge_genes_without_tags(const i_counter_t &cb_genes_count_by_id, ints_t &unmerged_cbs) const;

	names_t get_unmerged_names(const ints_t &unmerged_cbs) const;
	void parse_genes(const ints_t &unmerged_cbs, const s_counter_t &gene_counts, names_t &gene_names, ints_t &umis) const;
	i_counter_t count_cb_genes(bool logs = true) const;
	s_counter_t count_genes(const ints_t &unmerged_cbs) const;

	std::vector<double> get_reads_per_umis(const ints_t &unmerged_cbs) const;
	ints_t get_umig_coverage(const Estimator::i_counter_t &cb_genes_count_by_id) const;

	void print_text_output(const std::string &output_name, const names_t &gene_names,
						   const names_t &cell_names, const ints_t &umis) const;
	void print_bin_output(const std::string &bin_output_name, const IndropResult &results) const;

	IndropResult get_indrop_results(const CountMatrix cm, const i_counter_t &cb_genes_count_by_id,
											   const ints_t &unmerged_cbs, const ints_t &merge_n) const;

	std::string get_iseq_verbose(bam1_t *align_info, size_t read_prefix_length) const;
	std::string get_cb_top_verbose(const ints_t &unmerged_cbs) const;
	std::string get_cb_count_top_verbose(const i_counter_t &cb_genes_count_by_id) const;
	std::string get_genes_top_verbose(const s_counter_t &genes_count) const;

	void split_pairs(const SIHM &base, names_t &out1, ints_t out_2) const;
};
