#pragma once

#include "CellsDataContainer.h"
#include "Tools/ReadParameters.h"
#include "Tools/RefGenesContainer.h"

#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>

namespace Estimation
{
	class CellsDataContainer;

	class BamProcessor
	{
	private:
	public:
		BamProcessor(size_t read_prefix_length, const std::string &reads_params_names_str, const std::string &gtf_path);

	private:
		const bool _use_names_map;
		const size_t _read_prefix_length;
		const bool _use_genes_container;
		Tools::RefGenesContainer _genes_container;
		Tools::reads_params_map_t _reads_params;

	private:
		void parse_bam_file(const std::string &bam_name, bool print_result_bam, CellsDataContainer::s_i_map_t &cells_ids,
							CellsDataContainer::s_ii_hash_t &umig_cells_counts, CellsDataContainer &container) const;

		bool get_read_params(const std::string &read_name, Tools::ReadParameters &read_params) const;

		static int save_read_data(const std::string &chr_name, const std::string &cell_barcode, const std::string &umi,
							const std::string &gene, CellsDataContainer::s_i_map_t &cells_ids,
							CellsDataContainer::s_ii_hash_t &umig_cells_counts, CellsDataContainer &container);

		static std::string get_result_bam_name(const std::string &bam_name);

	public:
		CellsDataContainer::s_ii_hash_t parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
							 CellsDataContainer &container) const;

		std::string get_gene(const std::string &chr_name, BamTools::BamAlignment alignment) const;

		void write_alignment(BamTools::BamWriter &writer, BamTools::BamAlignment &alignment, const std::string &gene,
						 const Tools::ReadParameters &parameters) const;

		void fill_names_map(const std::string &reads_params_names_str);
	};
}