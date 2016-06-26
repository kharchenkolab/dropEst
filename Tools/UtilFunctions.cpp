#include "UtilFunctions.h"

#include <string.h>
#include <sstream>

namespace Tools
{
	unsigned edit_distance(const char *s1, const char *s2)
	{
		unsigned int s1len, s2len, x, y, lastdiag, olddiag;
		s1len = strlen(s1);
		s2len = strlen(s2);
		unsigned int column[s1len + 1];
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
				column[y] = MIN3(column[y] + 1, column[y - 1] + 1, lastdiag + (s1[y - 1] == s2[x - 1] ? 0 : 1));
				lastdiag = olddiag;
			}
		}
		return (column[s1len]);
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
}

