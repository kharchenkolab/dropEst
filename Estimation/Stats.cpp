#include "Stats.h"

namespace Estimation
{
	Stats::Stats()
	{}

	void Stats::add_merge_count(int count)
	{
		this->_merge_counts.push_back(count);
	}

	const Stats::int_list_t &Stats::get_merge_counts() const
	{
		return this->_merge_counts;
	}

	void Stats::inc(Stats::StringStatType counter, const std::string &name)
	{
		this->_named_counters[counter][name]++;
	}


	void Stats::get(Stats::StringStatType counter, Stats::str_list_t &names, Stats::int_list_t &counts) const
	{
		const s_cnt_t &s_counter = this->_named_counters[counter];
		for (s_cnt_t::const_iterator it = s_counter.begin(); it != s_counter.end(); ++it)
		{
			names.push_back(it->first);
			counts.push_back(it->second);
		}
	}

	Stats::int_list_t Stats::get(Stats::StringStatType counter) const
	{
		const s_cnt_t &s_counter = this->_named_counters[counter];
		int_list_t result;
		for (s_cnt_t::const_iterator it = s_counter.begin(); it != s_counter.end(); ++it)
		{
			result.push_back(it->second);
		}

		return result;
	}

	void Stats::inc_cells(CellStrStatType stat, const std::string &cell_barcode, const std::string &subtype)
	{
		this->_ss_cell_counters[stat][cell_barcode][subtype]++;
		this->_ss_cell_subtypes[stat].emplace(subtype);
	}

	void Stats::get_cells(CellStrStatType stat, str_list_t &types, str_list_t &subtypes, int_list_t &counts) const
	{
		std::copy(this->_ss_cell_subtypes[stat].begin(), this->_ss_cell_subtypes[stat].end(), std::back_inserter(subtypes));
		for (auto const &val : this->_ss_cell_counters[stat])
		{
			types.push_back(val.first);
			this->fill_by_types(val.second, subtypes, counts);
		}
	}

	void Stats::get_cells_filtered(CellStrStatType stat, const str_list_t &filter_barcodes, str_list_t &cell_barcodes,
	                               str_list_t &subtypes, int_list_t &counts) const
	{
		std::copy(this->_ss_cell_subtypes[stat].begin(), this->_ss_cell_subtypes[stat].end(), std::back_inserter(subtypes));
		for (auto const &barcode : filter_barcodes)
		{
			if (this->get_cell(stat, barcode, subtypes, counts))
			{
				cell_barcodes.push_back(barcode);
			}
		}
	}

	void Stats::fill_by_types(const s_cnt_t &counter, const str_list_t &types, int_list_t &counts) const
	{
		for (auto const &name : types)
		{
			s_cnt_t::const_iterator count = counter.find(name);
			counts.push_back(count == counter.end() ? 0 : count->second);
		}
	}

	void Stats::merge(const int_list_t &reassigned, const str_list_t &names)
	{
		for (int stat_num = 0; stat_num < CELL_S_STAT_SIZE; ++stat_num)
		{
			ss_cnt_t &cur_stat = this->_ss_cell_counters[stat_num];
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
				for (auto const &count: cell_from_it->second)
				{
					cell_to[count.first] += count.second;
				}

				cur_stat.erase(cell_from_it);
			}
		}
	}

	bool Stats::get_cell(CellStrStatType stat, const std::string &cell_barcode, const str_list_t &subtypes,
	                          int_list_t &counts) const
	{
		auto const &cur_stat = this->_ss_cell_counters[stat];
		ss_cnt_t::const_iterator counter_it = cur_stat.find(cell_barcode);
		if (counter_it == cur_stat.end())
			return false;

		this->fill_by_types(counter_it->second, subtypes, counts);
		return true;
	}

	const Stats::ss_cnt_t &Stats::get_raw_cell_stats(Stats::CellStrStatType stat) const
	{
		return this->_ss_cell_counters[stat];
	}
}