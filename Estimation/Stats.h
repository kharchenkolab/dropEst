#pragma once

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <boost/functional/hash.hpp>

class Stats
{
public:
	typedef std::vector<std::string> str_list_t;
	typedef std::vector<long> int_list_t;
	typedef std::vector<size_t> id_list_t;

	enum StringCounter
	{
		READS_BY_UMIG,
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

private:
	void get_chr_umi(const s_cnt_t &cell, const str_list_t &chr_names, int_list_t &counts) const;

public:
	Stats();

	void inc(StringCounter counter, const std::string &name);
	void get(StringCounter counter, str_list_t &names, int_list_t &counts) const;
	int_list_t get(StringCounter counter) const;

	void inc_cell_chr_umi(const std::string &chr_name, const std::string &cell_name, StatType type);
	void get_cell_chr_umi(StatType type, str_list_t &cell_names, str_list_t &chr_names, int_list_t &counts) const;
	void get_cell_chr_umi_filtered(StatType type, const str_list_t &filter_names, str_list_t &cell_names,
								   str_list_t &chr_names, int_list_t &counts) const;

	void add_merge_count(int count);
	const int_list_t& get_merge_counts() const;

	void merge(const int_list_t &reassigned, const str_list_t &names);
};