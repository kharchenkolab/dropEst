#pragma once

#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

#include "IndropResults.h"

namespace Estimation
{
	namespace Results
	{
		class IndropResultsWithoutUmi : public IndropResult
		{
			friend class boost::serialization::access;

		public:
			IndropResultsWithoutUmi()
			{ };

			IndropResultsWithoutUmi(const CountMatrix &cm, const Stats &stats, bool not_filtered)
				: IndropResult(cm, stats, std::vector<double>(), int_list_t(), not_filtered)
			{}

#ifdef R_LIBS
			virtual Rcpp::List get_main_r_vec(const std::string &filename) const override
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
									Named("merge.n") = wrap(this->merge_n),
									Named("exone_reads_by_cb") = wrap(this->exone_reads_by_cb),
									Named("fname") = wrap(filename));
			}
#endif
		};
	}
}