#include "Stats.h"

void Stats::inc_exone_chr_reads(std::string &chr_name)
{
	this->exone_chr_reads[chr_name]++;
}

void Stats::inc_nonexone_chr_reads(std::string &chr_name)
{
	this->nonexone_chr_reads[chr_name]++;
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
//	for (auto const &i : base)
	for (s_counter_t::const_iterator it = base.begin(); it != base.end(); ++it)
	{
		out_1.push_back(it->first);
		out_2.push_back(it->second);
	}
}

void Stats::add_merge_count(int count)
{
	this->merge_counts.push_back(count);
}

const Stats::int_list_t &Stats::get_merge_counts() const
{
	return this->merge_counts;
}

void Stats::inc_exone_cell_reads(std::string &cell_barcode)
{
	this->exone_cell_reads[cell_barcode]++;
}

void Stats::get_exone_cell_stats(Stats::str_list_t &exon_count_names, Stats::int_list_t &exon_counts) const
{
	Stats::split_pairs(this->exone_cell_reads, exon_count_names, exon_counts);
}

void Stats::inc_nonexone_cell_reads(std::string &cell_barcode)
{
	this->nonexone_cell_reads[cell_barcode]++;
}

void Stats::get_nonexone_cell_stats(Stats::str_list_t &nonexon_count_names, Stats::int_list_t &nonexon_counts) const
{
	Stats::split_pairs(this->nonexone_cell_reads, nonexon_count_names, nonexon_counts);
}
