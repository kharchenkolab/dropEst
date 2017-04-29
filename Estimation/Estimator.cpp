#include "Estimator.h"

#include "Estimation/BamProcessing/BamProcessor.h"
#include <Estimation/Results/IndropResults.h>
#include <Estimation/Results/BadCellsStats.h>
#include <Estimation/Results/IndropResultsWithoutUmi.h>
#include "Estimation/Merge/MergeStrategyFactory.h"
#include "Estimation/MergeUMIs/MergeUMIsStrategySimple.h"
#include "Tools/Logs.h"

#include <boost/range/adaptor/reversed.hpp>
#include <Estimation/BamProcessing/BamController.h>

using namespace std;

namespace Estimation
{
	const size_t Estimator::top_print_size;

	static bool comp_counters(const pair<string, int> &p1, const pair<string, int> &p2)
	{
		return p1.second > p2.second;
	}

	Estimator::Estimator(const boost::property_tree::ptree &config, bool merge_tags, const std::string &barcodes_filename)
		: merge_strategy(Merge::MergeStrategyFactory::get(config, merge_tags, barcodes_filename))
	{}

	Results::IndropResult Estimator::get_results(const CellsDataContainer &container, bool not_filtered, bool reads_output)
	{
		if (container.filtered_cells().empty())
		{
			L_WARN << "WARNING: filtered cells is empty. Maybe its too strict threshold or you forgot to run 'merge_and_filter'";
		}

		names_t gene_names(this->get_gene_names_sorted(container));

		Tools::trace_time("Compiling count matrix");
		names_t cell_names(this->get_filtered_cell_names(container));

		i_list_t count_matrix = this->get_count_matrix(gene_names, container, reads_output);

		Tools::trace_time("Done");

		Results::CountMatrix cm(cell_names, gene_names, count_matrix);
		return reads_output ? Results::IndropResultsWithoutUmi(cm, container, not_filtered)
			   : this->get_indrop_results(cm, container, not_filtered);
	}

	Results::BadCellsStats Estimator::get_bad_cells_results(const CellsDataContainer &container)
	{
		if (container.filtered_cells().empty())
		{
			L_WARN << "WARNING: filtered cells is empty. Maybe its too strict treshold or you forgot to run 'merge_and_filter'";
		}
		return Results::BadCellsStats();
	}

	Estimator::names_t Estimator::get_gene_names_sorted(const CellsDataContainer &genes_container) const
	{
		SIHM counter;
		for (size_t cell_index : genes_container.filtered_cells())
		{
			const CellsDataContainer::genes_t &cell_genes = genes_container.cell_genes(cell_index);
			for (CellsDataContainer::genes_t::const_iterator gm_it = cell_genes.begin(); gm_it != cell_genes.end(); ++gm_it)
			{
				counter[gm_it->first] += gm_it->second.size();
			}
		}

		L_TRACE << "\n" << counter.size() << " genes";

		s_counter_t gene_counts(counter.begin(), counter.end());
		sort(gene_counts.begin(), gene_counts.end(), comp_counters);

		L_TRACE << this->get_genes_top_verbose(gene_counts);

		names_t gene_names(gene_counts.size());
		for (int i = 0; i < gene_counts.size(); ++i)
		{
			gene_names[i] = gene_counts[i].first;
		}
		return gene_names;
	}

	Estimator::names_t Estimator::get_filtered_cell_names(const CellsDataContainer &genes_container) const
	{
		names_t cell_names;
		cell_names.reserve(genes_container.filtered_cells().size());
		for (auto const &cell_id : genes_container.filtered_cells())
		{
			cell_names.push_back(genes_container.cell_barcode(cell_id));
		}
		return cell_names;
	}

	Estimator::i_list_t Estimator::get_count_matrix(const names_t &gene_names, const CellsDataContainer &genes_container,
													bool reads_output) const
	{
		i_list_t count_matrix(genes_container.filtered_cells().size() * gene_names.size(), 0);
		std::unordered_map<std::string, size_t> gene_ids;

		for (size_t i = 0; i < gene_names.size(); ++i)
		{
			gene_ids[gene_names[i]] = i;
		}

		for (size_t col = 0; col < genes_container.filtered_cells().size(); col++)
		{
			for (auto const &gene_it : genes_container.cell_genes(genes_container.filtered_cells()[col]))
			{
				size_t row = gene_ids[gene_it.first];
				int cell_value = 0;

				if (reads_output)
				{
					for (auto const &umi : gene_it.second)
					{
						cell_value += umi.second.read_count;
					}
				}
				else
				{
					cell_value = gene_it.second.size();
				}

				count_matrix[(row * genes_container.filtered_cells().size()) + col] = cell_value;
			}
		}

		return count_matrix;
	}

	Results::IndropResult Estimator::get_indrop_results(const Results::CountMatrix &cm, const CellsDataContainer &genes_container,
														bool not_filtered) const
	{
		L_TRACE << "compiling diagnostic stats: ";

		doubles_t reads_per_umis(this->get_reads_per_umis(genes_container));
		L_TRACE << "reads/UMI";

		l_list_t umig_coverage(this->get_umig_coverage(genes_container));
		L_TRACE << "UMIg coverage";

		return Results::IndropResult(cm, genes_container, reads_per_umis, umig_coverage, not_filtered);
	}

	Estimator::doubles_t Estimator::get_reads_per_umis(const CellsDataContainer &genes_container) const
	{
		doubles_t reads_per_umis;
		reads_per_umis.reserve(genes_container.filtered_cells().size());

		for (size_t j = 0; j < genes_container.filtered_cells().size(); j++)
		{
			size_t umis_count = 0;
			double reads_per_umi = 0.0;
			for (auto const &gene_rec : genes_container.cell_genes(genes_container.filtered_cells()[j]))
			{
				for (auto const &umi_rec : gene_rec.second)
				{
					reads_per_umi += umi_rec.second.read_count;
				}

				umis_count += gene_rec.second.size();
			}
			reads_per_umi /= umis_count;
			reads_per_umis.push_back(reads_per_umi);
		}

		return reads_per_umis;
	}

	Estimator::l_list_t Estimator::get_umig_coverage(const CellsDataContainer &genes_container) const
	{
		l_list_t umig_coverage;
		s_set umigs_seen;
		for (const Tools::IndexedValue &gene_count : boost::adaptors::reverse(genes_container.cells_gene_counts_sorted()))
		{
			int new_umigs = 0;
			for (auto const &gene_rec : genes_container.cell_genes(gene_count.index))
			{
				for (auto const &umi_rec: gene_rec.second)
				{
					string umig = umi_rec.first + gene_rec.first;
					pair<s_set::const_iterator, bool> res = umigs_seen.emplace(umig);
					if (res.second)
					{
						new_umigs++;
					}
				}
			}
			umig_coverage.push_back(new_umigs);
		}
		return umig_coverage;
	}

	string Estimator::get_genes_top_verbose(const s_counter_t &genes) const
	{
		ostringstream ss;
		ss << "top genes:\n";
		for (size_t i = 0; i < min(genes.size(), Estimator::top_print_size); i++)
		{
			ss << genes[i].first << '\t' << genes[i].second << "\n";
		}
		return ss.str();
	}

	CellsDataContainer Estimator::get_cells_container(const names_t &files, bool bam_output, bool filled_bam,
													  const std::string &reads_params_names_str,
	                                                  const std::string &gtf_filename,
	                                                  CellsDataContainer::GeneMatchLevel gene_match_level)
	{
		CellsDataContainer container(this->merge_strategy, std::make_shared<MergeUMIs::MergeUMIsStrategySimple>(1), // TODO: Move 1 to parameter
		                             Estimator::top_print_size, gene_match_level);

		BamProcessing::BamController::parse_bam_files(files, bam_output, filled_bam, reads_params_names_str,
													  gtf_filename, container);

		container.set_initialized();

		container.merge_and_filter();
		return container;
	}
}
