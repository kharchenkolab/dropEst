#pragma once

#include "Estimation/CellsDataContainer.h"
#include "Tools/RefGenesContainer.h"

#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>

namespace Tools
{
	class ReadParameters;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class BamProcessor
		{
		protected:
			static const std::string GENE_TAG;
			static const std::string CB_TAG;
			static const std::string UMI_TAG;

		private:
			const size_t _read_prefix_length;
			Tools::RefGenesContainer _genes_container;


		private:
			void parse_bam_file(const std::string &bam_name, bool print_result_bam, CellsDataContainer::s_i_map_t &cells_ids,
								CellsDataContainer::s_uu_hash_t &umig_cells_counts, CellsDataContainer &container,
								long &total_reads, long &total_exonic_reads) const;

			void write_alignment(BamTools::BamWriter &writer, BamTools::BamAlignment &alignment, const std::string &gene,
								 const Tools::ReadParameters &parameters) const;

			std::string get_gene(const std::string &chr_name, BamTools::BamAlignment alignment) const;

			static size_t save_read_data(const std::string &chr_name, const std::string &cell_barcode, const std::string &umi,
									  const std::string &gene, CellsDataContainer::s_i_map_t &cells_ids,
									  CellsDataContainer::s_uu_hash_t &umig_cells_counts, CellsDataContainer &container);

			static std::string get_result_bam_name(const std::string &bam_name);

		protected:
			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) const;

		public:
			BamProcessor(size_t read_prefix_length, const std::string &gtf_path);

			CellsDataContainer::s_uu_hash_t parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
															CellsDataContainer &container) const;

			virtual void init_temporaries_before_parsing(bool save_read_name) const;

			virtual void release_temporaries_after_parsing() const;
		};
	}
}
