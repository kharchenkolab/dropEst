#pragma once

#include <string>

class TrimsCounter
{
public:
	enum StatType
	{
		NO_TRIM = 0,
		RC,
		POLY_A,
		A_TRIM,
		STAT_SIZE
	};

private:
	int stats[STAT_SIZE];
	std::string names[STAT_SIZE];

public:
	TrimsCounter();

	void inc(StatType type);
	int get(StatType type) const;

	std::string print()  const;
};
