#pragma once

#include <string>

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))

namespace Tools
{
	unsigned edit_distance(const char *s1, const char *s2);
	std::string reverse_complement(const std::string &s);
};
