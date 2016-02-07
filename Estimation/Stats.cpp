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

void Stats::get_cell_chr_umi(StatType type, str_list_t &cell_names, str_list_t &chr_names, int_list_t &counts) const
{
	std::copy(this->_chr_names.begin(), this->_chr_names.end(), std::back_inserter(chr_names));
	const ss_cnt_t &cur_stat = this->_cells_chr_umis_counts[type];
	for (ss_cnt_t::const_iterator cells_it = cur_stat.begin(); cells_it != cur_stat.end(); ++cells_it)
	{
		cell_names.push_back(cells_it->first);
		this->get_chr_umi(cells_it->second, chr_names, counts);
	}
}

void Stats::get_cell_chr_umi_filtered(StatType type, const str_list_t &filter_names, str_list_t &cell_names,
									  str_list_t &chr_names, int_list_t &counts) const
{
	std::copy(this->_chr_names.begin(), this->_chr_names.end(), std::back_inserter(chr_names));
	const ss_cnt_t &cur_stat = this->_cells_chr_umis_counts[type];
	for (str_list_t::const_iterator name_it = filter_names.begin(); name_it != filter_names.end(); ++name_it)
	{
		ss_cnt_t::const_iterator cell_it = cur_stat.find(*name_it);
		if (cell_it == cur_stat.end())
			continue;

		cell_names.push_back(*name_it);
		this->get_chr_umi(cell_it->second, chr_names, counts);
	}
}

void Stats::get_chr_umi(const s_cnt_t &cell, const str_list_t &chr_names, int_list_t &counts) const
{
	for (str_list_t::const_iterator name_it = chr_names.begin(); name_it != chr_names.end(); ++name_it)
	{
		s_cnt_t::const_iterator count = cell.find(*name_it);
		counts.push_back(count == cell.end() ? 0 : count->second);
	}
}

void Stats::merge(const int_list_t &reassigned, const str_list_t &names)
{
	ss_cnt_t &cur_stat = this->_cells_chr_umis_counts[EXONE];
	for (size_t ind = 0; ind < reassigned.size(); ++ind)
	{
		if (reassigned[ind] == ind)
			continue;

		const std::string &cur_name = names[ind];
		const std::string &target_name = names[reassigned[ind]];
		ss_cnt_t::const_iterator cell_from_it = cur_stat.find(cur_name);
		if (cell_from_it == cur_stat.end())
			continue;

		s_cnt_t &cell_to = cur_stat[target_name];
		for (s_cnt_t::const_iterator count_it = cell_from_it->second.begin();
			 count_it != cell_from_it->second.end(); ++count_it)
		{
			cell_to[count_it->first] += count_it->second;
		}

		cur_stat.erase(cell_from_it);
	}
}
