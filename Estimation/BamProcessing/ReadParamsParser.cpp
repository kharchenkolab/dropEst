#include <Tools/Logs.h>
#include <Tools/ReadParameters.h>
#include <Tools/RefGenesContainer.h>
#include "ReadParamsParser.h"
#include "BamController.h"

namespace Estimation
{
namespace BamProcessing
{
	ReadParamsParser::ReadParamsParser(const std::string &genes_filename, const BamTags &tags)
		: tags(tags)
	{
		if (genes_filename.length() != 0)
		{
			this->_genes_container = Tools::RefGenesContainer(genes_filename);
		}
	}

	bool ReadParamsParser::get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params)
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

	CellsDataContainer::Mark ReadParamsParser::get_gene(const std::string &chr_name, BamTools::BamAlignment alignment,
	                                                         std::string &gene) const
	{
		CellsDataContainer::Mark mark;
		gene = "";

		if (!this->_genes_container.is_empty())
			return this->get_gene_from_reference(chr_name, alignment, gene);

		if (!alignment.GetTag(this->tags.gene, gene))
		{
			gene = "";
			mark.add(CellsDataContainer::Mark::HAS_NOT_ANNOTATED);
		}
		else
		{
			mark.add(CellsDataContainer::Mark::HAS_EXONS);
		}

		return mark;
	}

	CellsDataContainer::Mark ReadParamsParser::get_gene_from_reference(const std::string &chr_name,
	                                                                    const BamTools::BamAlignment &alignment,
	                                                                    std::string &gene) const
	{
		CellsDataContainer::Mark mark;
		// TODO: parse CIGAR
		auto gene_set1 = _genes_container.get_gene_info(chr_name, alignment.Position, alignment.Position + 1);
		int end_position = alignment.GetEndPosition();
		auto gene_set2 = _genes_container.get_gene_info(chr_name, end_position - 1, end_position);

		if (gene_set1.empty() && gene_set2.empty())
			return mark;

		if (gene_set1.size() == 1 && gene_set2.size() == 1)
		{
			if (gene_set1.begin()->gene_name == gene_set2.begin()->gene_name)
			{
				mark.add(gene_set1.begin()->type);
				mark.add(gene_set2.begin()->type);

				gene = gene_set1.begin()->gene_name;
			}

			return mark;
		}

		if (gene_set1.size() <= 1 && gene_set2.size() <= 1)
		{
			auto const &non_empty = gene_set1.empty() ? *gene_set2.begin() : *gene_set1.begin();
			gene = non_empty.gene_name;
			mark.add(non_empty.type);
			mark.add(CellsDataContainer::Mark::HAS_NOT_ANNOTATED);

			return mark;
		}

		if (gene_set1.empty() || gene_set2.empty())
			return mark;

		Tools::RefGenesContainer::QueryResult exon_result1, exon_result2;
		if (!find_exon(gene_set1, exon_result1))
			return mark;

		if (!find_exon(gene_set2, exon_result2))
			return mark;

		if (!exon_result1.gene_name.empty() && !exon_result2.gene_name.empty())
		{
			if (exon_result1.gene_name != exon_result2.gene_name)
				return mark;

			gene = exon_result1.gene_name;

			mark.add(exon_result1.type);
			mark.add(exon_result2.type);
			return mark;
		}

		return mark;
	}

	bool ReadParamsParser::find_exon(Tools::RefGenesContainer::query_results_t query_results,
	                                  Tools::RefGenesContainer::QueryResult &exon_result) const
	{
		for (auto const &query_res : query_results)
		{
			if (query_res.type != Tools::GtfRecord::EXON)
				continue;

			if (exon_result.gene_name.empty())
			{
				exon_result = query_res;
				continue;
			}

			if (exon_result.gene_name != query_res.gene_name)
				return false;
		}

		return true;
	}

	bool ReadParamsParser::has_introns() const
	{
		return this->_genes_container.has_introns();
	}
}
}
