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
		typedef std::vector<size_t> ids_t;

		enum CellStatType
		{
			EXONE_READS_PER_CB,
			MERGES_COUNT_PER_CB,
			CELL_STAT_SIZE
		};

		enum CellStrStatType
		{
			EXONE_READS_PER_CHR_PER_CELL = 0,
			READS_PER_UMIG_PER_CELL,
			NON_EXONE_READS_PER_CHR_PER_CELL,
			READS_PER_UMI_PER_CELL,
			UMI_PER_CELL,
			CELL_S_STAT_SIZE
		};

		enum StrStrStatType
		{
			MERGE_INTERSECT_SIZE_BY_CELL,
			MERGE_REAL_INTERSECT_SIZE_BY_CELL,
			MERGE_EDIT_DISTANCE_BY_CELL,
			MERGE_REJECTION_BY_CELL,
			MERGE_PROB_BY_CELL,
			S_S_STAT_SIZE
		};

	public:
		typedef boost::unordered_map<std::string, int> s_cnt_t;
		typedef boost::unordered_set<std::string> str_set_t;
		typedef boost::unordered_map<std::string, s_cnt_t> ss_cnt_t;

	private:
		s_cnt_t _cell_counters[CELL_STAT_SIZE];

		ss_cnt_t _str_cell_counters[CELL_S_STAT_SIZE];
		str_set_t _str_cell_subtypes[CELL_S_STAT_SIZE];

		ss_cnt_t _str_str_counters[S_S_STAT_SIZE];

	private:
		void fill_by_types(const s_cnt_t &counter, const str_list_t &types, int_list_t &counts) const;

	public:
		Stats();

		void inc(CellStatType counter, const std::string &name);
		void get(CellStatType counter, str_list_t &names, int_list_t &counts) const;
		int_list_t get(CellStatType counter) const;
		const s_cnt_t& get_raw_stat(CellStatType stat) const;

		void add(StrStrStatType stat, const std::string &base_type, const std::string &subtype, long value);

		void inc(CellStrStatType stat, const std::string &cell_barcode, const std::string &subtype);
		void get(CellStrStatType stat, str_list_t &types, str_list_t &subtypes, int_list_t &counts) const;
		bool get(CellStrStatType stat, const std::string &cell_barcode, const str_list_t &subtypes,
				 int_list_t &counts) const;
		void get_filtered(CellStrStatType stat, const str_list_t &filter_barcodes, str_list_t &cell_barcodes,
						  str_list_t &subtypes, int_list_t &counts) const;

		const ss_cnt_t& get_raw(CellStrStatType stat) const;
		const ss_cnt_t& get_raw(StrStrStatType stat) const;

		void merge(const ids_t &reassigned, const str_list_t &cell_names);
	};
}