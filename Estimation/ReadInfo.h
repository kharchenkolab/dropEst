#pragma once

#include <Estimation/UMI.h>
#include <Tools/ReadParameters.h>
#include <string>

namespace Estimation
{
	class ReadInfo
	{
	public:
		static const char quality_offset = 64;

		const Tools::ReadParameters params;
		const std::string gene;
		const std::string chromosome_name;
		const UMI::Mark umi_mark;

		ReadInfo(const Tools::ReadParameters& params, const std::string &gene, const std::string &chromosome_name,
		         const UMI::Mark &umi_mark)
			: params(params)
			, gene(gene)
			, chromosome_name(chromosome_name)
			, umi_mark(umi_mark)
		{}
	};
}