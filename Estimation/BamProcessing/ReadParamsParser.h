#pragma once

#include <unordered_set>
#include <api/BamAlignment.h>

#include <Estimation/CellsDataContainer.h>
#include <Estimation/UMI.h>
#include <Tools/GeneAnnotation/RefGenesContainer.h>
#include <Tools/GeneAnnotation/RefGenesContainer.h>
#include "BamTags.h"

namespace Tools
{
	class ReadParameters;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadParamsParser
		{
		private:
			Tools::GeneAnnotation::RefGenesContainer _genes_container;
			bool _gene_in_chromosome_name;

		protected:
			const BamTags tags;

		private:
			static bool get_bam_tag(const BamTools::BamAlignment &alignment, const std::string &tag, std::string &value);

			UMI::Mark parse_read_type(const BamTools::BamAlignment &alignment, std::string &gene) const;
			UMI::Mark get_gene_from_reference(const std::string &chr_name, const BamTools::BamAlignment &alignment,
			                                  std::string &gene) const;

			bool find_exon(Tools::GeneAnnotation::RefGenesContainer::query_results_t query_results,
			               Tools::GeneAnnotation::RefGenesContainer::QueryResult &exon_result) const;

		public:
			ReadParamsParser(const std::string &genes_filename, const BamTags &tags, bool gene_in_chromosome_name);

			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params);
			UMI::Mark get_gene(const std::string &chr_name, BamTools::BamAlignment alignment,
			                                  std::string &gene) const;

			bool has_introns() const;
		};
	}
}
