#include "IndropResults.h"
#include "Stats.h"


IndropResult::IndropResult(const CountMatrix &cm, const std::vector<int> &non_exon_chr_counts,
						   const std::vector<std::string> &non_exon_chr_count_names, const std::vector<int> &exon_chr_counts,
						   const std::vector<std::string> &exon_chr_count_names, const std::vector<double> &reads_per_umi,
						   const std::vector<int> &umig_covered, const std::vector<int> &merge_n)
	: cm(cm)
	, non_exon_chr_counts(non_exon_chr_counts)
	, non_exon_chr_count_names(non_exon_chr_count_names)
	, reads_per_umi(reads_per_umi)
	, umig_covered(umig_covered)
	, exon_chr_counts(exon_chr_counts)
	, exon_chr_count_names(exon_chr_count_names)
	, merge_n(merge_n)
{}

IndropResult::IndropResult(const CountMatrix &cm, const Stats &stats, const std::vector<double> &reads_per_umi,
						   const std::vector<int> &umig_covered)
	: cm(cm)
	, reads_per_umi(reads_per_umi)
	, umig_covered(umig_covered)
	, merge_n(stats.get_merge_counts())
{
	stats.get_exone_cell_stats(this->exon_cell_count_tags, this->exon_cell_counts);
	stats.get_exone_chr_stats(this->exon_chr_count_names, this->exon_chr_counts);

	stats.get_nonexone_cell_stats(this->non_exon_cell_count_tags, this->non_exon_cell_counts);
	stats.get_nonexone_chr_stats(this->non_exon_chr_count_names, this->non_exon_chr_counts);
}

Rcpp::List IndropResult::get_r_table(const std::string &fname) const
{
	using namespace Rcpp;
	return List::create(Named("cell.names") = wrap(this->cm.cell_names),
						Named("gene.names") = wrap(this->cm.gene_names),
						Named("cm") = wrap(this->cm.counts),
						Named("none_chr_c") = wrap(this->non_exon_chr_counts),
						Named("none_chr_name") = wrap(this->non_exon_chr_count_names),
						Named("e_chr_c") = wrap(this->exon_chr_counts),
						Named("e_chr_name") = wrap(this->exon_chr_count_names),
						Named("none_cell_c") = wrap(this->non_exon_cell_counts),
						Named("none_cell_tag") = wrap(this->non_exon_cell_count_tags),
						Named("e_cell_c") = wrap(this->exon_cell_counts),
						Named("e_cell_tag") = wrap(this->exon_cell_count_tags),
						Named("rpu") = wrap(this->reads_per_umi),
						Named("umig.cov") = wrap(this->umig_covered),
						Named("merge.n") = wrap(this->merge_n),
						Named("fname") = wrap(fname));
}