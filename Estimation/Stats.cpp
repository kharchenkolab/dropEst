#include "Stats.h"

namespace Estimation
{
	Stats::id_set_t Stats::_presented_chromosomes[Stats::CHROMOSOME_STAT_SIZE];
	Stats::str_map_t Stats::_chromosome_inds;
	std::vector<std::string> Stats::_chromosome_names;

	Stats::Stats()
	{
		for (int i = 0; i < CellStatType::CELL_STAT_SIZE; ++i)
		{
			this->_stat_data[i] = 0;
		}
	}

	void Stats::inc(Stats::CellStatType type)
	{
		this->_stat_data[type]++;
	}

	void Stats::inc(CellChrStatType stat, const std::string &subtype)
	{
		auto chrom_it = Stats::_chromosome_inds.emplace(subtype, Stats::_chromosome_inds.size());
		if (chrom_it.second)
		{
			Stats::_chromosome_names.push_back(subtype);
		}

		Stats::_presented_chromosomes[stat].insert(chrom_it.first->second);
		this->_chromosome_stat_data[stat][chrom_it.first->second]++;
	}

	void Stats::merge(const Stats &source)
	{
		for (size_t i = 0; i < CellStatType::CELL_STAT_SIZE; ++i)
		{
			this->_stat_data[i] += source._stat_data[i];
		}

		for (size_t i = 0; i < CellChrStatType::CHROMOSOME_STAT_SIZE; ++i)
		{
			for (auto const &count: source._chromosome_stat_data[i])
			{
				this->_chromosome_stat_data[i][count.first] += count.second;
			}
		}
	}

	Stats::stat_t Stats::get(Stats::CellStatType type) const
	{
		return this->_stat_data[type];
	}

	bool Stats::get(CellChrStatType stat, stat_list_t &counts) const
	{
		auto const &cur_stat = this->_chromosome_stat_data[stat];
		if (cur_stat.empty())
			return false;

		for (auto const &id : Stats::_presented_chromosomes[stat])
		{
			auto count_it = cur_stat.find(id);
			counts.push_back(count_it == cur_stat.end() ? 0 : count_it->second);
		}

		return true;
	}

	Stats::str_list_t Stats::presented_chromosomes(Stats::CellChrStatType type)
	{
		str_list_t res;
		for (auto &chr_id : Stats::_presented_chromosomes[type])
		{
			res.push_back(Stats::_chromosome_names[chr_id]);
		}

		return res;
	}

	void Stats::dec(Stats::CellStatType type)
	{
		this->_stat_data[type]--;
	}
}
