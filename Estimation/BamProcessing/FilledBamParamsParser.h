#pragma once

#include "BamProcessor.h"
#include "ReadsParamsParser.h"

namespace Estimation
{
namespace BamProcessing
{
	class FilledBamParamsParser : public ReadsParamsParser
	{
	public:
		FilledBamParamsParser(const std::string &gtf_path);

		virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) override;
	};
}
}
