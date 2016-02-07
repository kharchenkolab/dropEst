#include "IndropResults.h"


IndropResult::IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
						   const std::vector<int> &umig_covered, bool not_filtered)
	: cm(cm)
	, reads_per_umi(reads_per_umi)
	, umig_covered(umig_covered)
	, merge_n(stats.get_merge_counts())
{
	if (not_filtered)
	{
		stats.get_cell_chr_umi(Stats::EXONE, this->ex_cell_names, this->chr_names, this->ex_cells_chr_umis_counts);
		this->chr_names.clear();
		stats.get_cell_chr_umi(Stats::NON_EXONE, this->nonex_cell_names, this->chr_names,
							   this->nonex_cells_chr_umis_counts);
	}
	else
	{
		stats.get_cell_chr_umi_filtered(Stats::EXONE, this->cm.cell_names, this->ex_cell_names, this->chr_names,
										this->ex_cells_chr_umis_counts);
		this->chr_names.clear();
		stats.get_cell_chr_umi_filtered(Stats::NON_EXONE, this->cm.cell_names, this->nonex_cell_names, this->chr_names,
							   this->nonex_cells_chr_umis_counts);
	}
}

#ifdef R_LIBS
Rcpp::List IndropResult::get_r_table(const std::string &fname) const
{
	using namespace Rcpp;

	return List::create(Named("cell.names") = wrap(this->cm.cell_names),
						Named("gene.names") = wrap(this->cm.gene_names),
						Named("cm") = wrap(this->cm.counts),
						Named("ex_cells_chr_counts") = wrap(this->ex_cells_chr_umis_counts),
						Named("nonex_cells_chr_counts") = wrap(this->nonex_cells_chr_umis_counts),
						Named("ex_counts_cell_names") = wrap(this->ex_cell_names),
						Named("nonex_counts_cell_names") = wrap(this->nonex_cell_names),
						Named("counts_chr_names") = wrap(this->chr_names),
						Named("rpu") = wrap(this->reads_per_umi),
						Named("umig.cov") = wrap(this->umig_covered),
						Named("merge.n") = wrap(this->merge_n),
						Named("fname") = wrap(fname));
}
#endif
