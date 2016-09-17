#pragma once

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_tree/ptree.hpp>

#include "CellsDataContainer.h"

namespace Estimation
{
	namespace Results
	{
		class CountMatrix;
		class IndropResult;
		class BadCellsStats;
	}

	class Estimator
	{
	private:
		typedef boost::unordered_map<std::string, int> SIHM;
		typedef boost::unordered_set<std::string> s_set;

		typedef std::vector<std::pair<std::string, int> > s_counter_t;
		typedef std::vector<double> doubles_t;
		typedef CellsDataContainer::names_t names_t;
		typedef std::vector<int> i_list_t;
		typedef std::vector<long> l_list_t;
		typedef CellsDataContainer::ids_t ids_t;

	public:
		typedef boost::unordered_map<std::string, boost::unordered_map<std::string, unsigned>> ss_u_hash_t; //Can't use long because of RInside (see BadCellsStats)

	private:
		static const size_t top_print_size = 10;

		const size_t read_prefix_length;
		const double min_merge_fraction;
		const unsigned max_merge_edit_distance;

		const unsigned min_genes_after_merge;
		const size_t barcode2_length;
		unsigned min_genes_before_merge;

		const std::string merge_type;

	public:
		Estimator(const boost::property_tree::ptree &config);

		Results::IndropResult get_results(const CellsDataContainer &container, bool not_filtered, bool reads_output);
		Results::BadCellsStats get_bad_cells_results(const CellsDataContainer &container);

		CellsDataContainer get_cells_container(const names_t &files, bool merge_tags, bool bam_output, bool filled_bam,
		                                       const std::string &reads_params_names_str, const std::string &gtf_filename,
		                                       const std::string &barcodes_filename);

	private:
		names_t get_filtered_cell_names(const CellsDataContainer &genes_container, const ids_t &unmerged_cells) const;

		i_list_t get_count_matrix(const ids_t &filtered_cells, const names_t &gene_names,
							   const CellsDataContainer &genes_container, bool reads_output) const;

		names_t get_gene_names_sorted(const CellsDataContainer &genes_container, const ids_t &filtered_cells) const;

		doubles_t get_reads_per_umis(const CellsDataContainer &genes_container, const ids_t &filtered_cells) const;

		l_list_t get_umig_coverage(const CellsDataContainer &genes_container) const;

		Results::IndropResult get_indrop_results(const Results::CountMatrix &cm, const CellsDataContainer &genes_conteiner,
										const ids_t &unmerged_cells, bool not_filtered) const;

		std::string get_genes_top_verbose(const s_counter_t &genes_count) const;
	};
}