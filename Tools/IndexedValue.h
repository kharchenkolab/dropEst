#pragma once

#include <cstddef>

namespace Tools
{
	struct IndexedValue
	{
		IndexedValue(size_t index, long value)
				: index(index), value(value)
		{}

		size_t index;
		long value;

		static bool value_less(const IndexedValue &ic1, const IndexedValue &ic2)
		{
			return ic1.value < ic2.value;
		}
	};
}