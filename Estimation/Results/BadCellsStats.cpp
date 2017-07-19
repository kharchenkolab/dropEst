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

			ss_cnt_t genes_reads, genes_umis;
			std::unordered_map<std::string, ss_cnt_t> reads_per_umig;
			for (auto cell_id : container.filtered_cells())
			{
				std::string cell_barcode = container.cell_barcode(cell_id);

				auto &cell_umis = genes_umis[cell_barcode];
				auto &cell_reads_p_umigs = reads_per_umig[cell_barcode];

				for (auto const &gene_umis : container.cell_genes(cell_id))
				{
					cell_umis[gene_umis.first] = (unsigned)gene_umis.second.size();
					auto &current_gene = genes_reads[cell_barcode][gene_umis.first];
					auto &out_gene_umis = cell_reads_p_umigs[gene_umis.first];

					for (auto const &umi_reads : gene_umis.second)
					{
						size_t read_count = umi_reads.second.read_count;
						current_gene += read_count;
						out_gene_umis [umi_reads.first] = (unsigned) read_count;
					}
				}
			}

			std::unordered_map<std::string, std::string> merge_targets;
			for (size_t cell_from_id = 0; cell_from_id < container.cell_barcodes_raw().size(); ++cell_from_id)
			{
				auto const &barcode_from = container.cell_barcodes_raw()[cell_from_id];
				auto const &barcode_to = container.cell_barcodes_raw()[container.merge_targets()[cell_from_id]];
				merge_targets[barcode_from] = barcode_to;
			}

			(*R)["d"] = List::create(
					Named("query_reads") = wrap(genes_reads),
		            Named("query_umis") = wrap(genes_umis),
		            Named("cell_gene_reads_per_chr") = wrap(container.stats().get_raw(Stats::GENE_READS_PER_CHR_PER_CELL)),
		            Named("cell_intergenic_reads_per_chr") = wrap(container.stats().get_raw(Stats::INTERGENIC_READS_PER_CHR_PER_CELL)),
		            Named("reads_per_umig") = wrap(reads_per_umig),
		            Named("umis_per_chr") = wrap(container.stats().get_raw(Stats::GENE_UMIS_PER_CHR_PER_CELL)),

					Named("has_exon_per_cell") = wrap(container.stats().get_raw(Stats::HAS_EXON_READS_PER_CB)),
					Named("has_intron_per_cell") = wrap(container.stats().get_raw(Stats::HAS_INTRON_READS_PER_CB)),
					Named("has_not_annotated_per_cell") = wrap(container.stats().get_raw(Stats::HAS_NOT_ANNOTATED_READS_PER_CB)),
					Named("total_reads_per_cell") = wrap(container.stats().get_raw(Stats::TOTAL_READS_PER_CB)),

		            Named("merges_count") = wrap(container.stats().get_raw(Stats::MERGES_COUNT_PER_CB)),

		            Named("merge_edit_distance") = wrap(container.stats().get_raw(Stats::MERGE_EDIT_DISTANCE_BY_CELL)),
		            Named("merge_probs") = wrap(container.stats().get_raw(Stats::MERGE_PROB_BY_CELL)),
		            Named("merge_targets") = merge_targets,
		            Named("excluded_cells") = wrap(container.excluded_cells())
			);

			R->parseEvalQ("saveRDS(d, '" + filename + "')");
			Tools::trace_time("Done");
		}
	}
}
