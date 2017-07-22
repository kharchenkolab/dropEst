#include "UtilFunctions.h"

namespace Tools
{
	unsigned edit_distance(const char *s1, const char *s2, bool skip_n)
	{
		size_t s1len, s2len;
		unsigned x, y, lastdiag, olddiag;
		s1len = strlen(s1);
		s2len = strlen(s2);
		unsigned column[s1len + 1];
		for (y = 1; y <= s1len; y++)
		{
			column[y] = y;
		}
		for (x = 1; x <= s2len; x++)
		{
			column[0] = x;
			for (y = 1, lastdiag = x - 1; y <= s1len; y++)
			{
				olddiag = column[y];
				int penalty = ((s1[y - 1] == s2[x - 1]) || (skip_n && (s1[y - 1] == 'N' || s2[x - 1] == 'N'))) ? 0 : 1;
				column[y] = MIN3(column[y] + 1, column[y - 1] + 1, lastdiag + penalty);
				lastdiag = olddiag;
			}
		}
		return column[s1len];
	}

	unsigned hamming_distance(const std::string &s1, const std::string &s2, bool skip_n)
	{
		if (s1.size() != s2.size())
			throw std::runtime_error("Strings should have equal length");

		unsigned ed = 0;
		for (std::string::size_type i = 0; i < s1.size(); ++i)
		{
			if (s1[i] != s2[i] && (!skip_n || (s1[i] != 'N' && s2[i] != 'N')))
			{
				++ed;
			}
		}

		return ed;
	}

	std::string reverse_complement(const std::string &s)
	{
		char rcs[s.length()];

		for (int i = 0; i < s.length(); i++)
		{
			switch (s[s.length() - i - 1])
			{
				case 'A':
					rcs[i] = 'T';
					break;
				case 'T':
					rcs[i] = 'A';
					break;
				case 'C':
					rcs[i] = 'G';
					break;
				case 'G':
					rcs[i] = 'C';
					break;
				default:
					rcs[i] = 'N';
					break;
			}
		}

		return std::string(rcs, s.length());
	}

	RInside* init_r()
	{
		RInside *r = RInside::instancePtr();
		if (r == nullptr)
			return new RInside(0, 0);

		return r;
	}
}

