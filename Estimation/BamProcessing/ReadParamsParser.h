#pragma once

#include <unordered_set>
#include <api/BamAlignment.h>

#include <Estimation/CellsDataContainer.h>
#include <Estimation/UMI.h>
#include <Tools/RefGenesContainer.h>
#include <Tools/RefGenesContainer.h>
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
			Tools::RefGenesContainer _genes_container;
			bool _gene_in_chromosome_name;
			std::unordered_set<std::string> _unexpected_read_types;

		protected:
			const BamTags tags;

		public:
			ReadParamsParser(const std::string &genes_filename, const BamTags &tags, bool gene_in_chromosome_name);

			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params);
			UMI::Mark get_gene(const std::string &chr_name, BamTools::BamAlignment alignment,
			                                  std::string &gene);

			bool find_exon(Tools::RefGenesContainer::query_results_t query_results,
			               Tools::RefGenesContainer::QueryResult &exon_result) const;

			bool has_introns() const;

			UMI::Mark get_gene_from_reference(const std::string &chr_name, const BamTools::BamAlignment &alignment,
			                                  std::string &gene) const;

			UMI::Mark parse_read_type(const BamTools::BamAlignment &alignment, std::string &gene) const;
		};
	}
}
