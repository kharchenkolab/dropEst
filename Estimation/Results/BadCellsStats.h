#pragma once

#include "Estimation/Estimator.h"

#include <RInside.h>
#include <boost/unordered_map.hpp>
#include <string>

namespace Estimation
{
	class Stats;
	namespace Results
	{
		class BadCellsStats
		{
		public:
			typedef Estimation::Estimator::ss_u_hash_t ss_cnt_t;

		private:
			ss_cnt_t genes_reads;
			ss_cnt_t genes_umis;

		public:

			void save_rds(const CellsDataContainer &container, const std::string &filename) const;

			BadCellsStats(const ss_cnt_t &genes_reads, const ss_cnt_t &genes_umis);
		};
	}
}
