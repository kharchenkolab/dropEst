#include "IndropResults.h"


IndropResult::IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
						   const std::vector<int> &umig_covered)
	: cm(cm)
	, reads_per_umi(reads_per_umi)
	, umig_covered(umig_covered)
	, merge_n(stats.get_merge_counts())
{
	stats.get_cell_chr_umi(this->cell_names, this->chr_names, this->cells_chr_umis_counts);
	stats.get_cell_chr_umi_exones_filtered(this->cm.cell_names, this->chr_names, this->filtered_cells_chr_umis_counts);
}

#ifdef R_LIBS
Rcpp::List IndropResult::get_r_table(const std::string &fname) const
{
	using namespace Rcpp;
	NumericVector vec_array(wrap(this->cells_chr_umis_counts));

	arma::cube cube_array(vec_array.begin(), Stats::ST_SIZE, this->cell_names.size(), this->chr_names.size());
	return List::create(Named("cell.names") = wrap(this->cm.cell_names),
						Named("gene.names") = wrap(this->cm.gene_names),
						Named("cm") = wrap(this->cm.counts),
						Named("cells_chr_counts") = wrap(cube_array),
						Named("filtered_cells_chr_counts") = wrap(this->filtered_cells_chr_umis_counts),
						Named("counts_cell_names") = wrap(this->cell_names),
						Named("counts_chr_names") = wrap(this->chr_names),
						Named("rpu") = wrap(this->reads_per_umi),
						Named("umig.cov") = wrap(this->umig_covered),
						Named("merge.n") = wrap(this->merge_n),
						Named("fname") = wrap(fname));
}
#endif
