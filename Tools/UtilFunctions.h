#pragma once

#include <string>
#include <vector>
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

	class PairHash
	{
	private:
		inline void hash_combine(std::size_t& seed) const
		{}

		template <typename T, typename... Rest>
		inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) const
		{
			std::hash<T> hasher;
			seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
			hash_combine(seed, rest...);
		}

	public:
		template <class T1, class T2>
		std::size_t operator () (const std::pair<T1,T2> &p) const
		{
			size_t h;
			this->hash_combine(h, p.first, p.second);
			return h;
		}
	};

	/// Estimate edit distance
	/// \param s1 string 1
	/// \param s2 string 2
	/// \param skip_n set all distances to N's equal to 0
	/// \param max_ed maximal edit distance used in the dynamic programming algorithm. Default: 10000.
	unsigned edit_distance(const char *s1, const char *s2, bool skip_n = true, unsigned max_ed=10000);
	unsigned hamming_distance(const std::string &s1, const std::string &s2, bool skip_n = true);
	double fpow(double base, long exp);
	RInside* init_r();

	std::string expand_tilde_in_path(const std::string &path);
	void copy_file(const std::string &src_file, const std::string &dest_file);
	std::string expand_relative_path(const std::string &src_fname, const std::string &rel_fname);
	std::string ltrim(const std::string &str);

	std::vector<size_t> parse_vec_from_string(const std::string &str, const std::string &split_pattern = " \t");
};
