#pragma once

#include <string>
#include <RInside.h>

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))

namespace Tools
{
	unsigned edit_distance(const char *s1, const char *s2, bool skip_n = true);
	std::string reverse_complement(const std::string &s);
	RInside* init_r();
};
