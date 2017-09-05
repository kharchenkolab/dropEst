#pragma once

#include "BamProcessor.h"
#include "ReadParamsParser.h"

namespace Estimation
{
namespace BamProcessing
{
	class FilledBamParamsParser : public ReadParamsParser
	{
	public:
		FilledBamParamsParser(const std::string &gtf_path, const BamTags &tags, bool gene_in_chromosome_name);

		virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) override;
	};
}
}
