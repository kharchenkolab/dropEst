#pragma once

#include <api/BamAlignment.h>
#include <Tools/RefGenesContainer.h>
#include <Estimation/CellsDataContainer.h>

namespace Tools
{
	class ReadParameters;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadsParamsParser
		{
		private:
			Tools::RefGenesContainer _genes_container;

		public:
			ReadsParamsParser(const std::string &genes_filename);

			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params);
			CellsDataContainer::Mark get_gene(const std::string &chr_name, BamTools::BamAlignment alignment,
			                                  std::string &gene) const;
		};
	}
}
