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
			public:
				virtual void merge(CellsDataContainer &container) const = 0;
			};
		}
	}
}