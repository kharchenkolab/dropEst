#pragma once

#include <string>
#include <vector>

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))

namespace Tools
{
	int edit_distance(const char *s1, const char *s2);

	std::vector<std::string> split(const std::string &s, char delim = ' ');
};