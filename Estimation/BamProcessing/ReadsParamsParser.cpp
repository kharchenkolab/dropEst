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

	CellsDataContainer::Mark ReadsParamsParser::get_gene(const std::string &chr_name, BamTools::BamAlignment alignment,
	                                                     std::string &gene) const
	{
		using Mark=CellsDataContainer::Mark;
		Mark mark;
		gene = "";

		if (!this->_genes_container.is_empty())
		{
			// TODO: parse CIGAR
			auto gene_set1 = this->_genes_container.get_gene_info(chr_name, alignment.Position, alignment.Position + 1);
			int end_position = alignment.GetEndPosition();
			auto gene_set2 = this->_genes_container.get_gene_info(chr_name, end_position - 1, end_position);

			if (gene_set1.size() > 1 || gene_set2.size() > 1)
				return mark;

			std::string gene1 = gene_set1.empty() ? "" : *gene_set1.begin();
			std::string gene2 = gene_set2.empty() ? "" : *gene_set2.begin();

			if (gene1.empty() && gene2.empty())
				return mark;

			mark.add(Mark::HAS_ANNOTATED);

			if (gene1.empty() || gene2.empty())
			{
				gene = gene1.empty() ? gene2 : gene1;
				mark.add(Mark::HAS_NOT_ANNOTATED);
				return mark;
			}

			if (gene1 != gene2)
				return mark;

			gene = gene1;
			return mark;
		}

		if (!alignment.GetTag(BamController::GENE_TAG, gene))
		{
			gene = "";
			return mark;
		}

		mark.add(Mark::HAS_ANNOTATED);
		return mark;
	}

	ReadsParamsParser::ReadsParamsParser(const std::string &genes_filename)
	{
		if (genes_filename.length() != 0)
		{
			this->_genes_container = Tools::RefGenesContainer(genes_filename);
		}
	}
}
}
