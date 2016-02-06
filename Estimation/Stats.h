#pragma once

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <boost/functional/hash.hpp>

class Stats
{
public:
	typedef std::vector<std::string> str_list_t;
	typedef std::vector<int> int_list_t;

	enum NamedCounter
	{
		NC_SIZE
	};

	enum StatType
	{
		EXONE = 0,
		NON_EXONE,
		ST_SIZE
	};

private:
	typedef boost::unordered_map<std::string, int> s_cnt_t;
	typedef boost::unordered_map<std::string, s_cnt_t> ss_cnt_t;
	typedef std::vector<ss_cnt_t> ss_typed_cnt_t;
	typedef boost::unordered_set<std::string> str_set_t;

	s_cnt_t _named_counters[NC_SIZE];

	int_list_t _merge_counts;

	ss_typed_cnt_t _cells_chr_umis_counts;
	str_set_t _chr_names;

public:
	Stats();

	void inc(NamedCounter counter, const std::string &name);
	void get(NamedCounter counter, str_list_t &names, int_list_t &counts) const;

	void inc_cell_chr_umi(const std::string &chr_name, const std::string &cell_name, StatType type);
	void get_cell_chr_umi(str_list_t &cell_names, str_list_t &chr_names, int_list_t &counts) const;

	void add_merge_count(int count);
	const int_list_t& get_merge_counts() const;
};

