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
		std::string cell_barcode;
		if (!alignment.GetTag(this->tags.cell_barcode, cell_barcode))
			return false;

		std::string umi;
		if (!alignment.GetTag(this->tags.umi, umi))
			return false;

		std::string cb_quality;
		alignment.GetTag(this->tags.cb_quality, cb_quality);

		std::string umi_quality;
		alignment.GetTag(this->tags.umi_quality, umi_quality);

		try
		{
			read_params = Tools::ReadParameters(cell_barcode, umi, cb_quality, umi_quality, this->_min_barcode_quality);
		}
		catch (std::runtime_error &error)
		{
			L_ERR << error.what();
			return false;
		}

		return true;
	}

	FilledBamParamsParser::FilledBamParamsParser(const std::string &gtf_path, const BamTags &tags,
	                                             bool gene_in_chromosome_name, int min_barcode_quality)
		: ReadParamsParser(gtf_path, tags, gene_in_chromosome_name)
		, _min_barcode_quality(min_barcode_quality)
	{}
}
}