#include "MultiSpacerOutcomesCounter.h"

#include <iomanip>
#include <sstream>

namespace TagsSearch
{
	void MultiSpacerOutcomesCounter::inc(StatType type)
	{
		++this->stats[type];
	}

	MultiSpacerOutcomesCounter::MultiSpacerOutcomesCounter(size_t spacer_num)
		: no_spacer(spacer_num, 0)
	{
		for (int i = 0; i < STAT_SIZE; ++i)
		{
			this->stats[i] = 0;
		}

		this->names[OK] = "OK";
		this->names[SHORT_SEQ] = "Short";
	}

	std::string MultiSpacerOutcomesCounter::print(double normalizer) const
	{
		std::ostringstream out_stream;
		out_stream << "Outcomes:\n[";
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << " (" << this->names[i] << ") ";
		}

		if (this->no_spacer.size() == 1)
		{
			out_stream << " (No spacer ) ";
		}
		else
		{
			for (size_t i = 0; i < this->no_spacer.size(); ++i)
			{
				out_stream << " (No spacer " << i << ") ";
			}
		}
		out_stream << "]" << std::endl;

		out_stream << "[";
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << this->stats[i] << " ";
		}
		for (int stat_val : this->no_spacer)
		{
			out_stream << stat_val << " ";
		}
		out_stream << "]" << std::endl;

		out_stream << "[" << std::setprecision(3);
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << 100.0 * this->stats[i] / normalizer << " ";
		}
		for (int stat_val : this->no_spacer)
		{
			out_stream << 100.0 * stat_val / normalizer << " ";
		}
		out_stream << "] %" << std::endl;

		return out_stream.str();
	}

	void MultiSpacerOutcomesCounter::inc_no_spacer(size_t spacer_ind)
	{
		this->no_spacer.at(spacer_ind)++;
	}
}