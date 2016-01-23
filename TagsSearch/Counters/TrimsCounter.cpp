#include "TrimsCounter.h"

#include <iomanip>
#include <sstream>

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
		stats[i] = 0;
	}

	names[NO_TRIM] = "no trim";
	names[RC] = "RC";
	names[POLY_A] = "PolyA";
	names[A_TRIM] = "-A";
}

std::string TrimsCounter::print() const
{
	std::ostringstream out_stream;
	out_stream << " trimst:[";
	for (int i = 0; i < STAT_SIZE; i++)
	{
		out_stream << " (" <<names[i] << ") ";
	}
	out_stream << "]" << std::endl;

	double normalizer = 0;
	out_stream << " trimst:[";
	for (int i = 0; i < STAT_SIZE; i++)
	{
		out_stream << stats[i] << " ";
		normalizer += stats[i];
	}
	out_stream << "]" << std::endl;

	out_stream << " trimst:[" << std::setprecision(3);
	for (int i = 0; i < STAT_SIZE; i++)
	{
		out_stream << 100.0 * stats[i] / normalizer << " ";
	}
	out_stream << "] %" << std::endl;

	return out_stream.str();
}
