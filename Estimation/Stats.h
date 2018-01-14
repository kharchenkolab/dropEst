#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Estimation
{
	class Stats
	{
	public:
		using str_list_t = std::vector<std::string>;
		using stat_list_t = std::vector<int>;
		using stat_t = int;

		enum CellStatType
		{
			TOTAL_READS_PER_CB,
			TOTAL_UMIS_PER_CB,
			CELL_STAT_SIZE
		};

		enum CellChrStatType
		{
			EXON_READS_PER_CHR_PER_CELL = 0,
			INTRON_READS_PER_CHR_PER_CELL,
			INTERGENIC_READS_PER_CHR_PER_CELL,
			CHROMOSOME_STAT_SIZE
		};

	private:
		using id_set_t = std::unordered_set<size_t>;
		using i_cnt_t = std::unordered_map<size_t, int>;
		using str_map_t = std::unordered_map<std::string, size_t>;
		using names_t = std::vector<std::string>;

	private:
		int _stat_data[CELL_STAT_SIZE];
		i_cnt_t _chromosome_stat_data[CHROMOSOME_STAT_SIZE];

		static id_set_t _presented_chromosomes[CHROMOSOME_STAT_SIZE];
		static str_map_t _chromosome_inds;
		static names_t _chromosome_names;

	public:
		Stats();

		void inc(CellStatType type);
		void dec(CellStatType type);
		void inc(CellChrStatType stat, const std::string &subtype);
		stat_t get(CellStatType type) const;
		bool get(CellChrStatType stat, stat_list_t &counts) const;

		void merge(const Stats &source);

		static size_t get_index(str_map_t &indexes, names_t &names, const std::string &type);
		static str_list_t presented_chromosomes(CellChrStatType type);
	};
}
