#include "MergeUMIsStrategyAbstract.h"

namespace Estimation
{
namespace Merge
{
namespace UMIs
{
	const std::string MergeUMIsStrategyAbstract::nucleotides = "ACGT";

	std::string MergeUMIsStrategyAbstract::fix_n_umi_with_random(const std::string &umi)
	{
		std::string target_umi(umi);
		for (char &c : target_umi)
		{
			if (c != 'N')
				continue;

			c = MergeUMIsStrategyAbstract::nucleotides[rand() % MergeUMIsStrategyAbstract::nucleotides.size()];
		}

		return target_umi;
	}
}
}
}