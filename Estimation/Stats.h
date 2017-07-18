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
			MERGES_COUNT_PER_CB,
			HAS_EXON_READS_PER_CB,
			HAS_INTRON_READS_PER_CB,
			HAS_NOT_ANNOTATED_READS_PER_CB,
			TOTAL_READS_PER_CB,
			CELL_STAT_SIZE
		};

		enum CellStrStatType
		{
			GENE_READS_PER_CHR_PER_CELL = 0,
			GENE_UMIS_PER_CHR_PER_CELL,
			INTERGENIC_READS_PER_CHR_PER_CELL,
			CELL_S_STAT_SIZE
		};

		enum StrStrStatType
		{
			MERGE_EDIT_DISTANCE_BY_CELL,
			S_S_STAT_SIZE
		};

		enum StrStrFloatType
		{
			MERGE_PROB_BY_CELL,
			S_S_FLOAT_SIZE
		};

	public:
		typedef boost::unordered_map<std::string, int> s_cnt_t;
		typedef boost::unordered_set<std::string> str_set_t;
		typedef boost::unordered_map<std::string, s_cnt_t> ss_cnt_t;
		typedef boost::unordered_map<std::string, boost::unordered_map<std::string, double>> ss_float_t;

	private:
		s_cnt_t _cell_counters[CELL_STAT_SIZE];

		ss_cnt_t _str_cell_counters[CELL_S_STAT_SIZE];
		str_set_t _str_cell_subtypes[CELL_S_STAT_SIZE];

		ss_cnt_t _str_str_counters[S_S_STAT_SIZE];
		ss_float_t _str_str_stats[S_S_FLOAT_SIZE];

	private:
		void fill_by_types(const s_cnt_t &counter, const str_list_t &types, int_list_t &counts) const;

	public:
		Stats();

		void inc(CellStatType counter, const std::string &name);
		void get(CellStatType counter, str_list_t &names, int_list_t &counts) const;
		int_list_t get(CellStatType counter) const;
		const s_cnt_t& get_raw(CellStatType stat) const;

		void set(StrStrStatType stat, const std::string &base_type, const std::string &subtype, int value);
		void set(StrStrFloatType stat, const std::string &base_type, const std::string &subtype, double value);

		void inc(CellStrStatType stat, const std::string &cell_barcode, const std::string &subtype);
		void dec(CellStrStatType stat, const std::string &cell_barcode, const std::string &subtype);
		void get(CellStrStatType stat, str_list_t &types, str_list_t &subtypes, int_list_t &counts) const;
		bool get(CellStrStatType stat, const std::string &cell_barcode, const str_list_t &subtypes,
				 int_list_t &counts) const;
		void get_filtered(CellStrStatType stat, const str_list_t &filter_barcodes, str_list_t &cell_barcodes,
						  str_list_t &subtypes, int_list_t &counts) const;

		const ss_cnt_t& get_raw(CellStrStatType stat) const;
		const ss_cnt_t& get_raw(StrStrStatType stat) const;
		const ss_float_t& get_raw(StrStrFloatType stat) const;

		void merge(const ids_t &reassigned, const str_list_t &cell_names);
	};
}