#pragma once

#include <string>

namespace TagsSearch
{
	class TwoBarcodesCounter
	{
	public:
		enum StatType
		{
			OK = 0,
			SHORT_READ1,
			SHORT_READ2,
			STAT_SIZE
		};

	private:
		int stats[STAT_SIZE];
		std::string names[STAT_SIZE];

	public:
		TwoBarcodesCounter();

		void inc(StatType type);
		int get(StatType type) const;

		std::string print(double normalizer) const;
	};
}