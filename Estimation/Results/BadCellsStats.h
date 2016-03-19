#pragma once

#include "IRTableProvider.h"
#include <Estimation/Stats.h>

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

		protected:
#ifdef R_LIBS
			void IndropResult::save_r_table(const std::string &filename) const override
			{
				using namespace Rcpp;

				RInside R(0, 0);
				R["d"] = List::create(Named("umig_cbs") = wrap(this->cm.cell_names),
				                    Named("umigs") = wrap(this->cm.gene_names),
				                    Named("umig_counts") = wrap(this->cm.counts),
				                    Named("gene_cbs") = wrap(this->ex_cells_chr_umis_counts),
				                    Named("genes") = wrap(this->nonex_cells_chr_umis_counts),
				                    Named("gene_counts") = wrap(this->ex_cell_names));
				L_TRACE << "writing R data to " << output_name;

				R.parseEvalQ(
						"d$umig_counts<-as.data.frame(matrix(d$umig_counts, length(d$umig_cbs), "
								"length(d$umigs), byrow = TRUE), row.names = d$umig_cbs); "
								"colnames(d$umig_counts)<-d$umigs; d$umigs<-NULL; d$umig_cbs<-NULL;");

				R.parseEvalQ(
						"d$gene_counts<-as.data.frame(matrix(d$gene_counts, length(d$gene_cbs), "
								"length(d$genes), byrow = TRUE), row.names = d$genes); "
								"colnames(d$gene_counts)<-d$genes; d$genes<-NULL; d$gene_cbs<-NULL;");


				R.parseEvalQ("saveRDS(d, '" + filename + "')");
				L_TRACE << "Done";
			}
#endif
		public:
			BadCellsStats(const Stats &stats)
			{
				stats.get_cells(Stats::UMIGS_READS_PER_CELL, this->umig_cell_barcodes, this->umigs, this->umigs_counts);
				stats.get_cells(Stats::GENES_READS_PER_CELL, this->genes_cell_barcodes, this->genes, this->genes_counts);
			}
		};
	}
}