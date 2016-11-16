#include "TwoBarcodesCounter.h"

#include <iomanip>
#include <sstream>

namespace TagsSearch
{
	void TwoBarcodesCounter::inc(StatType type)
	{
		++this->stats[type];
	}

	int TwoBarcodesCounter::get(StatType type) const
	{
		return this->stats[type];
	}

	TwoBarcodesCounter::TwoBarcodesCounter()
	{
		for (int i = 0; i < STAT_SIZE; ++i)
		{
			this->stats[i] = 0;
		}

		this->names[OK] = "OK";
		this->names[SHORT_READ1] = "Short read 1";
		this->names[SHORT_READ2] = "Short read 2";
	}

	std::string TwoBarcodesCounter::print(double normalizer) const
	{
		std::ostringstream out_stream;
		out_stream << "Outcomes:\n[";
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << " (" << this->names[i] << ") ";
		}
		out_stream << "]" << std::endl;

		out_stream << "[";
		for (int i = 0; i < STAT_SIZE; i++)
		{
			out_stream << this->stats[i] << " ";
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