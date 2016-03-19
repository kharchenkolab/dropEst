#pragma once

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <boost/functional/hash.hpp>

namespace Estimation
{
	class Stats
	{
	public:
		typedef std::vector<std::string> str_list_t;
		typedef std::vector<long> int_list_t;
		typedef std::vector<size_t> id_list_t;

		enum StringStatType
		{
			READS_BY_UMIG,
			READS_BY_CB,
			S_STAT_SIZE
		};

		enum CellStrStatType
		{
			EXONE_UMI_PER_CELL_PER_CHR = 0,
			NON_EXONE_UMI_PER_CELL_PER_CHR,
			UMIGS_READS_PER_CELL,
			GENES_READS_PER_CELL,
			CELL_S_STAT_SIZE
		};

	private:
		typedef boost::unordered_map<std::string, int> s_cnt_t;
		typedef boost::unordered_map<std::string, s_cnt_t> ss_cnt_t;
		typedef boost::unordered_set<std::string> str_set_t;

		s_cnt_t _named_counters[S_STAT_SIZE];

		int_list_t _merge_counts;

		ss_cnt_t _ss_cell_counters[CELL_S_STAT_SIZE];
		str_set_t _ss_cell_subtypes[CELL_S_STAT_SIZE];

	private:
		void fill_by_types(const s_cnt_t &counter, const str_list_t &types, int_list_t &counts) const;

	public:
		Stats();

		void inc(StringStatType counter, const std::string &name);
		void get(StringStatType counter, str_list_t &names, int_list_t &counts) const;
		int_list_t get(StringStatType counter) const;

		void inc_cells(CellStrStatType stat, const std::string &cell_barcode, const std::string &subtype);
		void get_cells(CellStrStatType stat, str_list_t &types, str_list_t &subtypes, int_list_t &counts) const;
		void get_cells_filtered(CellStrStatType stat, const str_list_t &filter_barcodes, str_list_t &cell_barcodes,
		                        str_list_t &subtypes, int_list_t &counts) const;

		void add_merge_count(int count);
		const int_list_t& get_merge_counts() const;

		void merge(const int_list_t &reassigned, const str_list_t &names);
	};
}