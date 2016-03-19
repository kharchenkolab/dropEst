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
			Stats::str_list_t umig_cell_barcodes;
			Stats::str_list_t umigs;
			Stats::int_list_t umigs_counts;
			Stats::str_list_t genes_cell_barcodes;
			Stats::str_list_t genes;
			Stats::int_list_t genes_counts;
			const Stats::ss_cnt_t &raw_umigs_stat;
			const Stats::ss_cnt_t &raw_genes_stat;

		protected:
#ifdef R_LIBS
			void save_r_table_old(const std::string &filename) const
			{
				using namespace Rcpp;

				L_TRACE << "writing R data to " << filename;
				RInside *R = RInside::instancePtr();
				if (R == nullptr)
				{
					R = new RInside(0, 0);
				}

				(*R)["d"] = List::create(Named("umig_cbs") = wrap(this->umig_cell_barcodes),
				                    Named("umigs") = wrap(this->umigs),
				                    Named("umig_counts") = wrap(this->umigs_counts),
				                    Named("gene_cbs") = wrap(this->genes_cell_barcodes),
				                    Named("genes") = wrap(this->genes),
				                    Named("gene_counts") = wrap(this->genes_counts));

				R->parseEvalQ(
						"d$umig_counts<-as.data.frame(matrix(d$umig_counts, length(d$umig_cbs), "
								"length(d$umigs), byrow = TRUE), row.names = d$umig_cbs); "
								"colnames(d$umig_counts)<-d$umigs; d$umigs<-NULL; d$umig_cbs<-NULL;");

				R->parseEvalQ(
						"d$gene_counts<-as.data.frame(matrix(d$gene_counts, length(d$gene_cbs), "
								"length(d$genes), byrow = TRUE), row.names = d$genes); "
								"colnames(d$gene_counts)<-d$genes; d$genes<-NULL; d$gene_cbs<-NULL;");


				R->parseEvalQ("saveRDS(d, '" + filename + "')");
				L_TRACE << "Done";
			}

			virtual void save_r_table(const std::string &filename) const override
			{
				using namespace Rcpp;

				L_TRACE << "writing R data to " << filename;
				RInside *R = RInside::instancePtr();
				if (R == nullptr)
				{
					R = new RInside(0, 0);
				}

				(*R)["umigs"] = this->raw_umigs_stat;
				(*R)["genes"] = this->raw_genes_stat;

				R->parseEval("for(n in names(umigs)) {umigs[[n]] <- as.list(umigs[[n]])}");
				R->parseEval("for(n in names(genes)) {genes[[n]] <- as.list(genes[[n]])}");
				R->parseEval("d <- list(umigs=umigs, genes=genes)");

				R->parseEvalQ("saveRDS(d, '" + filename + "')");
				L_TRACE << "Done";
			}
#endif
		public:
			BadCellsStats(const Stats &stats)
				: raw_genes_stat(stats.get_raw_cell_stats(Stats::GENES_READS_PER_CELL))
				, raw_umigs_stat(stats.get_raw_cell_stats(Stats::UMIGS_READS_PER_CELL))
			{
//				stats.get_cells(Stats::UMIGS_READS_PER_CELL, this->umig_cell_barcodes, this->umigs, this->umigs_counts);
//				stats.get_cells(Stats::GENES_READS_PER_CELL, this->genes_cell_barcodes, this->genes, this->genes_counts);
			}
		};
	}
}
