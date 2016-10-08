#include <Tools/Logs.h>
#include <Tools/ReadParameters.h>
#include "ReadsParamsParser.h"
#include "BamController.h"

namespace Estimation
{
namespace BamProcessing
{
	bool ReadsParamsParser::get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params)
	{
		try
		{
			read_params = Tools::ReadParameters(alignment.Name);
		}
		catch (std::runtime_error &error)
		{
			L_ERR << error.what();
			return false;
		}

		return true;
	}

	std::string ReadsParamsParser::get_gene(const std::string &chr_name, BamTools::BamAlignment alignment) const
	{
		if (!this->_genes_container.is_empty())
			return this->_genes_container.get_gene_info(chr_name, alignment.Position,
														alignment.Position + alignment.Length).name();

		std::string gene;
		if (!alignment.GetTag(BamController::GENE_TAG, gene))
			return "";

		return gene;
	}

	ReadsParamsParser::ReadsParamsParser(const std::string &gtf_path)
	{
		if (gtf_path.length() != 0)
		{
			this->_genes_container = Tools::RefGenesContainer(gtf_path);
		}
	}
}
}