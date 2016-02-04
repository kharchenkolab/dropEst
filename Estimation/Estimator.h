#pragma once

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_tree/ptree.hpp>

#include "GenesContainer.h"
#include "IndropResults.h"
#include "Stats.h"

class Estimator
{
private:
	typedef boost::unordered_map<std::string, int> SIHM;
	typedef boost::unordered_set<std::string> s_set;

	typedef std::vector<std::pair<std::string, int> > s_counter_t;
	typedef std::vector<double> doubles_t;
	typedef GenesContainer::names_t names_t;
	typedef GenesContainer::ints_t ints_t;
	typedef GenesContainer::ids_t ids_t;

private:
	static const size_t top_print_size;

	size_t read_prefix_length;
	double min_merge_fraction;
	int max_merge_edit_distance;

	int min_genes_after_merge;
	int min_genes_before_merge;

public:
	Estimator(const boost::property_tree::ptree &config);

	IndropResult get_results(const names_t &files, bool merge_tags);

private:
	names_t get_unmerged_names(const GenesContainer &genes_container, const ids_t &unmerged_cells) const;
	void fill_table(const ids_t &unmerged_cells, const s_counter_t &gene_counts,
					const GenesContainer &genes_container, names_t &gene_names_header, ints_t &umis_table) const;

	s_counter_t count_genes(const GenesContainer &genes_container, const ids_t &unmerged_cells) const;

	doubles_t get_reads_per_umis(const GenesContainer &genes_container, const ids_t &unmerged_cells) const;
	ints_t get_umig_coverage(const GenesContainer &genes_container) const;

	IndropResult get_indrop_results(const CountMatrix cm, const GenesContainer &genes_conteiner,
									const ids_t &unmerged_cells) const;

	std::string get_cb_top_verbose(const GenesContainer &genes_container, const ids_t &unmerged_cells) const;
	std::string get_genes_top_verbose(const s_counter_t &genes_count) const;
};
