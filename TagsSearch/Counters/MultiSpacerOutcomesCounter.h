#pragma once

#include <string>
#include <vector>

namespace TagsSearch
{
	class MultiSpacerOutcomesCounter
	{
	public:
		enum StatType
		{
			OK = 0,
			SHORT_SEQ,
			STAT_SIZE
		};

	private:
		int stats[STAT_SIZE];
		std::string names[STAT_SIZE];
		std::vector<int> no_spacer;

	public:
		MultiSpacerOutcomesCounter(size_t spacer_num);

		void inc(StatType type);
		void inc_no_spacer(size_t spacer_ind);

		std::string print(double normalizer)  const;
	};
}