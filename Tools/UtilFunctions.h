#pragma once

#include <string>
#include <RInside.h>

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))

namespace Tools
{
	class ReverseComplement
	{
	private:
		char complements[std::numeric_limits<char>::max()];
	public:
		ReverseComplement();
		std::string rc(const std::string &s) const;
	};
	unsigned edit_distance(const char *s1, const char *s2, bool skip_n = true, unsigned max_ed=10000);
	unsigned hamming_distance(const std::string &s1, const std::string &s2, bool skip_n = true);
	double fpow(double base, long exp);
	RInside* init_r();
};
