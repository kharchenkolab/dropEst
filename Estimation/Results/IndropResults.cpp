#include "IndropResults.h"

#include <RInside.h>

#include <Estimation/CellsDataContainer.h>
#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>

#include "boost/filesystem.hpp"

namespace Estimation
{
	namespace Results
	{
		IndropResult::IndropResult(const CountMatrix &cm, const CellsDataContainer &container,
								   const std::vector<double> &reads_per_umi, const int_list_t &umig_covered, bool not_filtered)
				: cm(cm)
				, reads_per_umi(reads_per_umi)
				, umig_covered(umig_covered)
		{
			for (size_t cell_ind = 0; cell_ind < container.cell_barcodes_raw().size(); ++cell_ind)
			{
				if (container.is_cell_merged(cell_ind))
					continue;

				const auto &cb = container.cell_barcode(cell_ind);
				size_t sum_reads = 0;
				for (auto const& gene : container.cell_genes(cell_ind))
				{
					for (auto const &umi : gene.second)
					{
						this->reads_by_umig_cbs.push_back(cb);
						this->reads_by_umig_umis.push_back(umi.first);
						this->reads_by_umig.push_back(umi.second);
						sum_reads += umi.second;
					}
				}

				this->exone_reads_by_cb.push_back(sum_reads);
			}

			if (not_filtered)
			{
				L_TRACE << "Fill exone results";
				container.stats().get(Stats::EXONE_READS_PER_CHR_PER_CELL, this->ex_cell_names, this->ex_chr_names,
						  this->ex_cells_chr_reads_counts);

				L_TRACE << "Fill nonexone results";
				container.stats().get(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, this->nonex_cell_names, this->nonex_chr_names,
						  this->nonex_cells_chr_reads_counts);
			}
			else
			{
				L_TRACE << "Fill exone results";
				container.stats().get_filtered(Stats::EXONE_READS_PER_CHR_PER_CELL, this->cm.cell_names, this->ex_cell_names,
								   this->ex_chr_names, this->ex_cells_chr_reads_counts);

				L_TRACE << "Fill nonexone results";
				container.stats().get_filtered(Stats::NON_EXONE_READS_PER_CHR_PER_CELL, this->cm.cell_names,
								   this->nonex_cell_names, this->nonex_chr_names, this->nonex_cells_chr_reads_counts);
			}
		}

		void IndropResult::save_rds(const CellsDataContainer &container, const std::string &filename) const
		{
			using namespace Rcpp;

			RInside *R = Tools::init_r();
			(*R)["data"] = this->get_main_r_vec(filename);
			Tools::trace_time("Writing R data to " + filename);

			R->parseEvalQ(
					"data$ex_cells_chr_counts<-as.data.frame(matrix(data$ex_cells_chr_counts, length(data$ex_counts_cell_names), "
							"length(data$ex_chr_names), byrow = TRUE), row.names = data$ex_counts_cell_names); "
							"colnames(data$ex_cells_chr_counts)<-data$ex_chr_names; data$ex_counts_cell_names<-NULL; data$ex_chr_names<-NULL");

			R->parseEvalQ(
					"data$nonex_cells_chr_counts<-as.data.frame(matrix(data$nonex_cells_chr_counts, length(data$nonex_counts_cell_names), "
							"length(data$nonex_chr_names), byrow = TRUE), row.names = data$nonex_counts_cell_names); "
							"colnames(data$nonex_cells_chr_counts)<-data$nonex_chr_names; data$nonex_counts_cell_names<-NULL;"
							"data$nonex_chr_names<-NULL;");

			R->parseEvalQ("data$cm <- matrix(data$cm, length(data$gene.names), length(data$cell.names), byrow=T);"
								  "rownames(data$cm) <- data$gene.names; colnames(data$cm) <- data$cell.names");

			(*R)["report_data"] = this->get_report_r_vec(container);
//			R->parseEval((std::string)"source('" + PROJ_BIN_PATH + "/Builder.R', echo=TRUE, verbose=TRUE)");
			R->parseEval((std::string)"source('" + PROJ_BIN_PATH + "/Builder.R')");
//			boost::filesystem::remove("Report.log");
//			boost::filesystem::remove("Report.aux");
//			boost::filesystem::remove("Report.tex");
//			boost::filesystem::remove_all("figures");

			R->parseEvalQ("saveRDS(data, '" + filename + "')");
			Tools::trace_time("Done");
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
						 Named("ex_chr_names") = wrap(this->ex_chr_names),
						 Named("nonex_chr_names") = wrap(this->nonex_chr_names),
						 Named("rpu") = wrap(this->reads_per_umi),
						 Named("umig.cov") = wrap(this->umig_covered),
						 Named("reads_by_umig") = wrap(this->reads_by_umig),
						 Named("reads_by_umig_cbs") = wrap(this->reads_by_umig_cbs),
						 Named("reads_by_umig_umis") = wrap(this->reads_by_umig_umis),
						 Named("exone_reads_by_cb") = wrap(this->exone_reads_by_cb),
						 Named("fname") = wrap(filename));
		}

		Rcpp::List IndropResult::get_report_r_vec(const CellsDataContainer &container) const
		{
			using namespace Rcpp;
			return List::create(Named("merge_type") = wrap(container.merge_type()),
								Named("merge_probs") = wrap(container.stats().get_raw(Stats::MERGE_PROB_BY_CELL)),
								Named("report_script") = wrap(PROJ_BIN_PATH + (std::string)"/Report.R")
			);
		}
	}
}
