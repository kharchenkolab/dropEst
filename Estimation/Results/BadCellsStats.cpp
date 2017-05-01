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

				auto &cell_umis = genes_umis[cell_barcode];
				auto &cell_reads_p_umis = reads_per_umi[cell_barcode];
				auto &cell_reads_p_umigs = reads_per_umig[cell_barcode];

				for (auto const &gene_umis : container.cell_genes(cell_id))
				{
					cell_umis[gene_umis.first] = (unsigned)gene_umis.second.size();
					auto &current_gene = genes_reads[cell_barcode][gene_umis.first];
					for (auto const &umi_reads : gene_umis.second)
					{

						size_t read_count = umi_reads.second.read_count;
						current_gene += read_count;
						cell_reads_p_umis[umi_reads.first] += read_count;
						cell_reads_p_umigs[umi_reads.first + gene_umis.first] = (unsigned) read_count;
					}
				}
			}

			(*R)["d"] = List::create(
					Named("genes_reads") = wrap(genes_reads),
		            Named("genes_umis") = wrap(genes_umis),
		            Named("cell_exone_reads_per_chr") = wrap(container.stats().get_raw(Stats::GENE_READS_PER_CHR_PER_CELL)),
		            Named("cell_nonexone_reads_per_chr") = wrap(container.stats().get_raw(Stats::INTERGENIC_READS_PER_CHR_PER_CELL)),
		            Named("reads_per_umig") = wrap(reads_per_umig),
		            Named("umis_per_chr") = wrap(container.stats().get_raw(Stats::GENE_UMIS_PER_CHR_PER_CELL)),

					Named("has_exon_per_cell") = wrap(container.stats().get_raw(Stats::HAS_EXON_READS_PER_CB)),
					Named("has_intron_per_cell") = wrap(container.stats().get_raw(Stats::HAS_INTRON_READS_PER_CB)),
					Named("has_not_annotated_per_cell") = wrap(container.stats().get_raw(Stats::HAS_NOT_ANNOTATED_READS_PER_CB)),

		            Named("merges_count") = wrap(container.stats().get_raw(Stats::MERGES_COUNT_PER_CB)),

		            Named("merge_edit_distance") = wrap(container.stats().get_raw(Stats::MERGE_EDIT_DISTANCE_BY_CELL)),
		            Named("merge_probs") = wrap(container.stats().get_raw(Stats::MERGE_PROB_BY_CELL)),
		            Named("merge_targets") = wrap(container.stats().get_raw(Stats::MERGE_TARGET_BY_BASE)),
		            Named("excluded_cells") = wrap(container.excluded_cells()),
		            Named("umi_per_cell") = wrap(reads_per_umi)
			);

			R->parseEvalQ("saveRDS(d, '" + filename + "')");
			Tools::trace_time("Done");
		}
	}
}
