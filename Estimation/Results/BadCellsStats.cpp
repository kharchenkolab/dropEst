#include "BadCellsStats.h"

#include "Estimation/Stats.h"
#include "Tools/Logs.h"
#include <Rcpp.h>
#include <Tools/UtilFunctions.h>

namespace Estimation
{
	namespace Results
	{
		void BadCellsStats::save_rds(const CellsDataContainer &container, const std::string &filename) const
		{
			using namespace Rcpp;

			Tools::trace_time("Writing R data to " + filename);
			RInside *R = Tools::init_r();

			ss_cnt_t genes_reads, genes_umis, reads_per_umig, reads_per_umi;
			for (auto cell_id : container.filtered_cells())
			{
				std::string cell_barcode = container.cell_barcode(cell_id);
				auto &cell = genes_umis[cell_barcode];
				auto &cell_umis = reads_per_umi[cell_barcode];
				auto &cell_umigs = reads_per_umig[cell_barcode];
				size_t sum_reads = 0;

				for (auto const &gene : container.cell_genes(cell_id))
				{
					cell[gene.first] = (unsigned)gene.second.size();
					auto &current_gene = genes_reads[cell_barcode][gene.first];
					for (auto const &umi : gene.second)
					{
						current_gene += umi.second;
						cell_umis[umi.first] += umi.second;
						cell_umigs[umi.first + gene.first] = (unsigned)umi.second;
					}
				}
			}

			(*R)["d"] = List::create(
					Named("genes_reads") = wrap(genes_reads),
		            Named("genes_umis") = wrap(genes_umis),
		            Named("cell_exone_reads_per_chr") = wrap(container.stats().get_raw(Stats::EXONE_READS_PER_CHR_PER_CELL)),
		            Named("cell_nonexone_reads_per_chr") = wrap(container.stats().get_raw(Stats::NON_EXONE_READS_PER_CHR_PER_CELL)),
		            Named("reads_per_umig") = wrap(reads_per_umig),
		            Named("umis_per_chr") = wrap(container.stats().get_raw(Stats::EXONE_UMIS_PER_CHR_PER_CELL)),

		            Named("merges_count") = wrap(container.stats().get_raw_stat(Stats::MERGES_COUNT_PER_CB)),

		            Named("merge_edit_distance") = wrap(container.stats().get_raw(Stats::MERGE_EDIT_DISTANCE_BY_CELL)),
		            Named("merge_rejections") = wrap(container.stats().get_raw(Stats::MERGE_REJECTION_BY_CELL)),
		            Named("merge_probs") = wrap(container.stats().get_raw(Stats::MERGE_PROB_BY_CELL)),
		            Named("excluded_cells") = wrap(container.excluded_cells()),
		            Named("umi_per_cell") = wrap(reads_per_umi)
			);

			R->parseEvalQ("saveRDS(d, '" + filename + "')");
			Tools::trace_time("Done");
		}
	}
}
