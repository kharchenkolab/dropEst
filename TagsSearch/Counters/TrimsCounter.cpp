#include "TrimsCounter.h"

#include <iomanip>
#include <sstream>

namespace TagsSearch
{
	void TrimsCounter::inc(StatType type)
	{
		++this->stats[type];
	}

	int TrimsCounter::get(StatType type) const
	{
		return this->stats[type];
	}

	TrimsCounter::TrimsCounter()
	{
		for (int i = 0; i < STAT_SIZE; ++i)
		{
			this->stats[i] = 0;
		}

		this->names[NO_TRIM] = "no trim";
		this->names[RC] = "RC";
		this->names[POLY_A] = "PolyA";
		this->names[A_TRIM] = "-A";
	}

	std::string TrimsCounter::print() const
	{
		std::ostringstream out_stream;
		out_stream << "\tTrim stat:\n[";
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << " (" << this->names[i] << ") ";
		}
		out_stream << "]" << std::endl;

		double normalizer = 0;
		out_stream << "[";
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << this->stats[i] << " ";
			normalizer += this->stats[i];
		}
		out_stream << "]" << std::endl;

		out_stream << "[" << std::setprecision(3);
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << 100.0 * this->stats[i] / normalizer << " ";
		}
		out_stream << "] %" << std::endl;

		return out_stream.str();
	}
}