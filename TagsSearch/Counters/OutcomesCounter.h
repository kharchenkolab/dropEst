#pragma once

#include <string>

namespace TagsSearch
{
	class OutcomesCounter
	{
	public:
		enum StatType
		{
			OK = 0,
			NO_SPACER,
			SHORT_SEQ,
			SPACER_MODIFIED,
			STAT_SIZE
		};

	private:
		int stats[STAT_SIZE];
		std::string names[STAT_SIZE];

	public:
		OutcomesCounter();

		void inc(StatType type);
		int get(StatType type) const;

		std::string print(double normalizer)  const;
	};
}