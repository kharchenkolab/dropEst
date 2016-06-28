#pragma once

#include "IRTableProvider.h"
#include <Estimation/Stats.h>

#ifdef R_LIBS
#include <RInside.h>
#endif

namespace Estimation
{
	namespace Results
	{
		class BadCellsStats : public IRTableProvider
		{	
		private:
			Stats::ss_cnt_t genes_umis;
			Stats::ss_cnt_t genes_reads;
			Stats::ss_cnt_t cell_exone_reads_per_chr;
			Stats::ss_cnt_t cell_nonexone_reads_per_chr;
			Stats::ss_cnt_t reads_per_umig;
			Stats::ss_cnt_t reads_per_umi;

			Stats::ss_cnt_t merge_edit_distance;
			Stats::ss_cnt_t merge_intersect_size;
			Stats::ss_cnt_t merge_rejections;

			Stats::s_cnt_t merges_count;
			std::vector<std::string> excluded_cells;

		protected:
#ifdef R_LIBS
			virtual void save_r_table(const std::string &filename) const override;
#endif
		public:
			BadCellsStats(const Stats &stats, const Stats::ss_cnt_t &genes_reads, const Stats::ss_cnt_t &genes_umis,
			              const std::vector<std::string> &excluded_cells);
		};
	}
}
