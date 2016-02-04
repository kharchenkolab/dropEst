#pragma once

#include <string>
#include <boost/unordered_map.hpp>

#include <boost/functional/hash.hpp>

class Stats
{
public:
	typedef std::vector<std::string> str_list_t;
	typedef std::vector<int> int_list_t;

private:
	typedef boost::unordered_map<std::string, int> s_counter_t;

	s_counter_t exone_chr_reads;
	s_counter_t nonexone_chr_reads;

	s_counter_t exone_cell_reads;
	s_counter_t nonexone_cell_reads;

	int_list_t merge_counts;


private:
	static void split_pairs(const s_counter_t &base, str_list_t &out_1, int_list_t &out_2);

public:
	void inc_exone_chr_reads(std::string &chr_name);
	void get_exone_chr_stats(str_list_t &exon_count_names, int_list_t &exon_counts) const;

	void inc_nonexone_chr_reads(std::string &chr_name);
	void get_nonexone_chr_stats(str_list_t &nonexon_count_names, int_list_t &nonexon_counts) const;

	void inc_exone_cell_reads(std::string &cell_barcode);
	void get_exone_cell_stats(str_list_t &exon_count_names, int_list_t &exon_counts) const;

	void inc_nonexone_cell_reads(std::string &cell_barcode);
	void get_nonexone_cell_stats(str_list_t &nonexon_count_names, int_list_t &nonexon_counts) const;

	void add_merge_count(int count);
	const int_list_t& get_merge_counts() const;
};

