#include "FilledBamParamsParser.h"
#include "BamController.h"

#include <Tools/Logs.h>
#include <Tools/ReadParameters.h>
#include <api/BamAlignment.h>

namespace Estimation
{
namespace BamProcessing
{
	bool FilledBamParamsParser::get_read_params(const BamTools::BamAlignment &alignment,
											 Tools::ReadParameters &read_params)
	{
		std::string barcode;
		if (!alignment.GetTag(BamController::CB_TAG, barcode))
			return false;

		std::string umi;
		if (!alignment.GetTag(BamController::UMI_TAG, umi))
			return false;

		try
		{
			read_params = Tools::ReadParameters(alignment.Name, barcode, umi);
		}
		catch (std::runtime_error &error)
		{
			L_ERR << error.what();
			return false;
		}

		return true;
	}

	FilledBamParamsParser::FilledBamParamsParser(const std::string &gtf_path, int gene_match_level)
			: ReadsParamsParser(gtf_path, gene_match_level)
	{}
}
}