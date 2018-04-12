#pragma once

#include <Estimation/CellsDataContainer.h>

namespace Estimation
{
	namespace Merge
	{
		namespace UMIs
		{
			class MergeUMIsStrategyAbstract
			{
			protected:
				static const std::string nucleotides;

			protected:
				static std::string fix_n_umi_with_random(const std::string &umi);
			public:
				virtual void merge(CellsDataContainer &container) const = 0;
			};
		}
	}
}