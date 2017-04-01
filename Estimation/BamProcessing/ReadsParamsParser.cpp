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
		{
			// TODO: parse CIGAR
			std::string gene1 = this->_genes_container.get_gene_info(chr_name, alignment.Position, alignment.Position + 1).name();
			int end_position = alignment.GetEndPosition();
			std::string gene2 = this->_genes_container.get_gene_info(chr_name, end_position - 1, end_position).name();

			if (this->_gene_match_level == GeneMatchLevel::BOTH)
				return (gene1 == gene2) ? gene1 : "";

			if (this->_gene_match_level == GeneMatchLevel::ONE)
				return (gene1 == gene2) ? "" : (gene1 == "" ? gene2 : "");

			if (gene1 == "")
				return gene2;

			if (gene2 == "")
				return gene1;

			if (gene1 != gene2)
				return "";

			return gene1;
		}

		std::string gene;
		if (!alignment.GetTag(BamController::GENE_TAG, gene))
			return "";

		return gene;
	}

	ReadsParamsParser::ReadsParamsParser(const std::string &genes_filename, int gene_match_level)
		: _gene_match_level(gene_match_level)
	{
		if (genes_filename.length() != 0)
		{
			this->_genes_container = Tools::RefGenesContainer(genes_filename);
		}
	}
}
}
