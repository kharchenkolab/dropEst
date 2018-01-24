#pragma once

#include "BamProcessor.h"
#include "ReadParamsParser.h"

namespace Estimation
{
namespace BamProcessing
{
	class FilledBamParamsParser : public ReadParamsParser
	{
	protected:
		const int _min_barcode_quality;

	public:
		FilledBamParamsParser(const std::string &gtf_path, const BamTags &tags, bool gene_in_chromosome_name,
		                      int min_barcode_quality);
		bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) override;
	};
}
}
