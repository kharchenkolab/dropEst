#include "IndropResults.h"

#include <Tools/Logs.h>
#include <RInside.h>

namespace Estimation
{
	namespace Results
	{
		IndropResult::IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
		                           const int_list_t &umig_covered, bool not_filtered)
				: cm(cm), reads_per_umi(reads_per_umi), umig_covered(umig_covered), merge_n(stats.get_merge_counts()),
				  reads_by_umig(stats.get(Stats::READS_BY_UMIG)), exone_reads_by_cb(stats.get(Stats::EXONE_READS_PER_CB))
		{
			if (not_filtered)
			{
				L_TRACE << "Fill exone results";
				stats.get_cells(Stats::EXONE_READS_PER_CHR_PER_CELL, this->ex_cell_names, this->chr_names,
				                this->ex_cells_chr_reads_counts);
				this->chr_names.clear();
				L_TRACE << "Fill nonexone results";
				stats.get_cells(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, this->nonex_cell_names, this->chr_names,
				                this->nonex_cells_chr_reads_counts);
			}
			else
			{
				L_TRACE << "Fill exone results";
				stats.get_cells_filtered(Stats::EXONE_READS_PER_CHR_PER_CELL, this->cm.cell_names, this->ex_cell_names,
				                         this->chr_names, this->ex_cells_chr_reads_counts);
				this->chr_names.clear();
				L_TRACE << "Fill nonexone results";
				stats.get_cells_filtered(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, this->cm.cell_names,
				                         this->nonex_cell_names, this->chr_names, this->nonex_cells_chr_reads_counts);
			}
		}

		void IndropResult::save_rds(const std::string &filename) const
		{
			using namespace Rcpp;

			RInside *R = RInside::instancePtr();
			if (R == nullptr)
			{
				R = new RInside(0, 0);
			}
			(*R)["d"] = this->get_main_r_vec(filename);
			L_TRACE << "writing R data to " << filename;

			R->parseEvalQ(
					"d$ex_cells_chr_counts<-as.data.frame(matrix(d$ex_cells_chr_counts, length(d$ex_counts_cell_names), "
							"length(d$counts_chr_names), byrow = TRUE), row.names = d$ex_counts_cell_names); "
							"colnames(d$ex_cells_chr_counts)<-d$counts_chr_names; d$ex_counts_cell_names<-NULL;");

			R->parseEvalQ(
					"d$nonex_cells_chr_counts<-as.data.frame(matrix(d$nonex_cells_chr_counts, length(d$nonex_counts_cell_names), "
							"length(d$counts_chr_names), byrow = TRUE), row.names = d$nonex_counts_cell_names); "
							"colnames(d$nonex_cells_chr_counts)<-d$counts_chr_names; d$nonex_counts_cell_names<-NULL;"
							"d$counts_chr_names<-NULL;");

			R->parseEvalQ("saveRDS(d, '" + filename + "')");
			L_TRACE << "Done";
		}

		Rcpp::List IndropResult::get_main_r_vec(const std::string &filename) const
		{
			using namespace Rcpp;
			return List::create(Named("cell.names") = wrap(this->cm.cell_names),
						 Named("gene.names") = wrap(this->cm.gene_names),
						 Named("cm") = wrap(this->cm.counts),
						 Named("ex_cells_chr_counts") = wrap(this->ex_cells_chr_reads_counts),
						 Named("nonex_cells_chr_counts") = wrap(this->nonex_cells_chr_reads_counts),
						 Named("ex_counts_cell_names") = wrap(this->ex_cell_names),
						 Named("nonex_counts_cell_names") = wrap(this->nonex_cell_names),
						 Named("counts_chr_names") = wrap(this->chr_names),
						 Named("rpu") = wrap(this->reads_per_umi),
						 Named("umig.cov") = wrap(this->umig_covered),
						 Named("merge.n") = wrap(this->merge_n),
						 Named("reads_by_umig") = wrap(this->reads_by_umig),
						 Named("exone_reads_by_cb") = wrap(this->exone_reads_by_cb),
						 Named("fname") = wrap(filename));
		}
	}
}
