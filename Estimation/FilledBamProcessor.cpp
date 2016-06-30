#include <Tools/Logs.h>
#include "FilledBamProcessor.h"

namespace Estimation
{
	bool FilledBamProcessor::get_read_params(const BamTools::BamAlignment &alignment,
	                                         Tools::ReadParameters &read_params) const
	{
		std::string barcode;
		if (!alignment.GetTag(BamProcessor::CB_TAG, barcode))
			return false;

		std::string umi;
		if (!alignment.GetTag(BamProcessor::UMI_TAG, umi))
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

	FilledBamProcessor::FilledBamProcessor(size_t read_prefix_length, const std::string &gtf_path)
		: BamProcessor(read_prefix_length, gtf_path)
	{}


}