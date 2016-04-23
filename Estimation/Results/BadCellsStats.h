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
			const Stats::ss_cnt_t &genes_umis;
			const Stats::ss_cnt_t &genes_reads;
			const Stats::ss_cnt_t &cell_exone_reads_per_chr;
			const Stats::ss_cnt_t &cell_nonexone_reads_per_chr;

			const Stats::s_cnt_t &merges_count;
			const Stats::s_cnt_t &reads_per_cb;

		protected:
#ifdef R_LIBS
			virtual void save_r_table(const std::string &filename) const override;
#endif
		public:
			BadCellsStats(const Stats &stats, const Stats::ss_cnt_t &genes_reads, const Stats::ss_cnt_t &genes_umis);
		};
	}
}
