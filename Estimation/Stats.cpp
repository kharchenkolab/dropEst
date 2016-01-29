#include "Stats.h"

void Stats::inc_exone_chr_reads(std::string &chr_name)
{
	this->exone_chr_reads[chr_name]++;
//	if (!this->exone_chr_reads.count(chr_name))
//	{
//		this->exone_chr_reads[chr_name] = 1;
//	}
//	else
//	{
//		this->exone_chr_reads[chr_name]++;
//	}
}

void Stats::inc_nonexone_chr_reads(std::string &chr_name)
{
	this->nonexone_chr_reads[chr_name]++;
//	if (!this->nonexone_chr_reads.count(chr_name))
//	{
//		this->nonexone_chr_reads[chr_name] = 1;
//	}
//	else
//	{
//		this->nonexone_chr_reads[chr_name]++;
//	}
}

void Stats::get_exone_chr_stats(Stats::str_list_t &exon_count_names, Stats::int_list_t &exon_counts) const
{
	Stats::split_pairs(this->exone_chr_reads, exon_count_names, exon_counts);
}

void Stats::get_nonexone_chr_stats(Stats::str_list_t &nonexon_count_names, Stats::int_list_t &nonexon_counts) const
{
	Stats::split_pairs(this->nonexone_chr_reads, nonexon_count_names, nonexon_counts);
}

void Stats::split_pairs(const Stats::s_counter_t &base, Stats::str_list_t &out_1, Stats::int_list_t &out_2)
{
	for (auto const &i : base)
	{
		out_1.push_back(i.first);
		out_2.push_back(i.second);
	}
}

void Stats::add_merge_count(long count)
{
	this->merge_counts.push_back(count);
}

const Stats::int_list_t &Stats::get_merge_counts() const
{
	return this->merge_counts;
}
