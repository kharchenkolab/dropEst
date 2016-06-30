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

		std::string gene;
		if (!alignment.GetTag(BamProcessor::GENE_TAG, gene))
		{
			gene = "";
		}

//	read_params.ReadParameters()
		throw std::runtime_error("Not implemented");
		return true;
	}

	FilledBamProcessor::FilledBamProcessor(size_t read_prefix_length, const std::string &gtf_path)
		: BamProcessor(read_prefix_length, gtf_path)
	{}


}