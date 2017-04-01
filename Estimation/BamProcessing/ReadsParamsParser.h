#pragma once

#include <api/BamAlignment.h>
#include <Tools/RefGenesContainer.h>

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
			bool _exons_only;

		public:
			ReadsParamsParser(const std::string &genes_filename, bool exons_only = false);

			virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params);
			std::string get_gene(const std::string &chr_name, BamTools::BamAlignment alignment) const;
		};
	}
}
