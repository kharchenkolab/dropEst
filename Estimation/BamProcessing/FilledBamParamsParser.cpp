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
		if (!alignment.GetTag(this->tags.cb, barcode))
			return false;

		std::string umi;
		if (!alignment.GetTag(this->tags.umi, umi))
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

	FilledBamParamsParser::FilledBamParamsParser(const std::string &gtf_path, const BamTags &tags)
		: ReadParamsParser(gtf_path, tags)
	{}
}
}