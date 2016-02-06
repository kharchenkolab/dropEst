#include "Stats.h"

Stats::Stats()
{
	this->_cells_chr_umis_counts.resize(ST_SIZE);
}

void Stats::add_merge_count(int count)
{
	this->_merge_counts.push_back(count);
}

const Stats::int_list_t &Stats::get_merge_counts() const
{
	return this->_merge_counts;
}

void Stats::inc(Stats::NamedCounter counter, const std::string &name)
{
	this->_named_counters[counter][name]++;
}


void Stats::get(Stats::NamedCounter counter, Stats::str_list_t &names, Stats::int_list_t &counts) const
{
	const s_cnt_t &s_counter = this->_named_counters[counter];
	for (s_cnt_t::const_iterator it = s_counter.begin(); it != s_counter.end(); ++it)
	{
		names.push_back(it->first);
		counts.push_back(it->second);
	}
}

void Stats::inc_cell_chr_umi(const std::string &chr_name, const std::string &cell_name, StatType type)
{
	this->_cells_chr_umis_counts[type][cell_name][chr_name]++;
	this->_chr_names.emplace(chr_name);
}

void Stats::get_cell_chr_umi(Stats::str_list_t &cell_names, Stats::str_list_t &chr_names,
							 Stats::int_list_t &counts) const
{
	std::copy(this->_chr_names.begin(), this->_chr_names.end(), chr_names.begin());
	for (int stat_type = 0; stat_type < ST_SIZE; ++stat_type)
	{
		const ss_cnt_t &cur_stat = this->_cells_chr_umis_counts[stat_type];
		for (ss_cnt_t::const_iterator cells_it = cur_stat.begin(); cells_it != cur_stat.end(); ++cells_it)
		{
			cell_names.push_back(cells_it->first);
			for (str_list_t::const_iterator name_it = chr_names.begin(); name_it != chr_names.end(); ++name_it)
			{
				s_cnt_t::const_iterator count = cells_it->second.find(*name_it);
				counts.push_back(count == cells_it->second.end() ? 0 : count->second);
			}
		}
	}
}
