#include "Estimator.h"

#include "Estimation/BamProcessing/BamProcessor.h"
#include "Estimation/BamProcessing/BamProcessorFactory.h"
#include <Estimation/Results/IndropResults.h>
#include <Estimation/Results/BadCellsStats.h>
#include <Estimation/Results/IndropResultsWithoutUmi.h>
#include "Estimation/Merge/MergeStrategyFactory.h"
#include "Tools/Logs.h"

#include <boost/range/adaptor/reversed.hpp>

using namespace std;

namespace Estimation
{
	const size_t Estimator::top_print_size;

	static bool comp_counters(const pair<string, int> &p1, const pair<string, int> &p2)
	{
		return p1.second > p2.second;
	}

	Estimator::Estimator(const boost::property_tree::ptree &config)
		: read_prefix_length(config.get<size_t>("read_prefix_length"))
		, min_merge_fraction(config.get<double>("min_merge_fraction"))
		, max_merge_edit_distance(config.get<unsigned>("max_merge_edit_distance"))
		, min_genes_after_merge(config.get<unsigned>("min_genes_after_merge"))
		, barcode2_length(config.get<size_t>("barcode2_length", 0))
		, min_genes_before_merge(config.get<unsigned>("min_genes_before_merge"))
	{
		if (this->min_genes_after_merge > 0 && this->min_genes_after_merge < this->min_genes_before_merge)
		{
			this->min_genes_before_merge = this->min_genes_after_merge;
		}

		if (barcode2_length == 0)
		{
			L_WARN << "Barcode2 length is equal to 0";
		}
	}

	Results::IndropResult Estimator::get_results(const CellsDataContainer &container, bool not_filtered, bool reads_output)
	{
		ids_t filtered_cells = container.filtered_cells();
		if (filtered_cells.empty())
		{
			L_WARN << "WARNING: filtered cells is empty. Maybe its too strict threshold or you forgot to run 'merge_and_filter'";
		}

		L_TRACE << this->get_cb_top_verbose(container, filtered_cells);

		L_TRACE << filtered_cells.size() << " valid (with >=" << this->min_genes_after_merge << " genes) cells with ";

		names_t gene_names(this->get_gene_names_sorted(container, filtered_cells));

		Tools::trace_time("Compiling count matrix");
		names_t cell_names(this->get_filtered_cell_names(container, filtered_cells));

		ints_t count_matrix = this->get_count_matrix(filtered_cells, gene_names, container, reads_output);

		Tools::trace_time("Done");

		Results::CountMatrix cm(cell_names, gene_names, count_matrix);
		return reads_output ? Results::IndropResultsWithoutUmi(cm, container.stats(), not_filtered)
			   : this->get_indrop_results(cm, container, filtered_cells, not_filtered);
	}

	Results::BadCellsStats Estimator::get_bad_cells_results(const CellsDataContainer &container)
	{
		if (container.filtered_cells().empty())
		{
			L_WARN << "WARNING: filtered cells is empty. Maybe its too strict treshold or you forgot to run 'merge_and_filter'";
		}
		return Results::BadCellsStats(this->get_reads_per_genes_per_cells_count(container),
									  this->get_umis_per_genes_per_cells_count(container), container.excluded_cells());
	}

	Estimator::names_t Estimator::get_gene_names_sorted(const CellsDataContainer &genes_container,
														const ids_t &unmerged_cells) const
	{
		SIHM counter;
		for (size_t cell_index: unmerged_cells)
		{
			const CellsDataContainer::genes_t &cell_genes = genes_container.cell_genes(cell_index);
			for (CellsDataContainer::genes_t::const_iterator gm_it = cell_genes.begin(); gm_it != cell_genes.end(); ++gm_it)
			{
				counter[gm_it->first] += gm_it->second.size();
			}
		}

		L_TRACE << counter.size() << " genes";

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

	Estimator::names_t Estimator::get_filtered_cell_names(const CellsDataContainer &genes_container,
														  const ids_t &unmerged_cells) const
	{
		names_t cell_names;
		cell_names.reserve(unmerged_cells.size());
		for (auto const &cell_id : unmerged_cells)
		{
			cell_names.push_back(genes_container.cell_barcode(cell_id));
		}
		return cell_names;
	}

	Estimator::ints_t Estimator::get_count_matrix(const ids_t &unmerged_cells, const names_t &gene_names,
									  const CellsDataContainer &genes_container, bool reads_output) const
	{
		ints_t count_matrix(unmerged_cells.size() * gene_names.size(), 0);
		std::unordered_map<std::string, size_t> gene_ids;

		for (size_t i = 0; i < gene_names.size(); ++i)
		{
			gene_ids[gene_names[i]] = i;
		}

		for (size_t col = 0; col < unmerged_cells.size(); col++)
		{
			for (auto const &gene_it : genes_container.cell_genes(unmerged_cells[col]))
			{
				size_t row = gene_ids[gene_it.first];
				long cell_value = 0;

				if (reads_output)
				{
					for (auto const &umi : gene_it.second)
					{
						cell_value += umi.second;
					}
				}
				else
				{
					cell_value = gene_it.second.size();
				}

				count_matrix[(row * unmerged_cells.size()) + col] = cell_value;
			}
		}

		return count_matrix;
	}

	Results::IndropResult Estimator::get_indrop_results(const Results::CountMatrix &cm, const CellsDataContainer &genes_container,
											   const ids_t &unmerged_cells, bool not_filtered) const
	{
		L_TRACE << "compiling diagnostic stats: ";

		doubles_t reads_per_umis(this->get_reads_per_umis(genes_container, unmerged_cells));
		L_TRACE << "reads/UMI";

		ints_t umig_coverage(this->get_umig_coverage(genes_container));
		L_TRACE << "UMIg coverage";

		return Results::IndropResult(cm, genes_container.stats(), reads_per_umis, umig_coverage, not_filtered);
	}

	Estimator::doubles_t Estimator::get_reads_per_umis(const CellsDataContainer &genes_container,
													   const ids_t &unmerged_cells) const
	{
		doubles_t reads_per_umis;
		reads_per_umis.reserve(unmerged_cells.size());

		for (size_t j = 0; j < unmerged_cells.size(); j++)
		{
			size_t umis_count = 0;
			double reads_per_umi = 0.0;
			for (auto const &gene_rec : genes_container.cell_genes(unmerged_cells[j]))
			{
				for (auto const &umi_rec : gene_rec.second)
				{
					reads_per_umi += umi_rec.second;
				}

				umis_count += gene_rec.second.size();
			}
			reads_per_umi /= umis_count;
			reads_per_umis.push_back(reads_per_umi);
		}

		return reads_per_umis;
	}

	Estimator::ints_t Estimator::get_umig_coverage(const CellsDataContainer &genes_container) const
	{
		ints_t umig_coverage;
		s_set umigs_seen;
		for (const Tools::IndexedValue &gene_count : boost::adaptors::reverse(genes_container.cells_genes_counts_sorted()))
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


	string Estimator::get_cb_top_verbose(const CellsDataContainer &genes_container, const ids_t &unmerged_cells) const
	{
		stringstream ss;
		if (unmerged_cells.size() > 0)
		{
			ss << "top CBs:\n";
			for (size_t i = 0; i < min(unmerged_cells.size(), Estimator::top_print_size); i++)
			{
				ss << genes_container.cell_genes(unmerged_cells[i]).size() << "\t" <<
						genes_container.cell_barcode(unmerged_cells[i]) << "\n";
			}
		}
		else
		{
			ss << "no valid CBs found\n";
		}

		return ss.str();
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

	CellsDataContainer Estimator::get_cells_container(const names_t &files, bool merge_tags, bool bam_output,
	                                                  bool filled_bam, const std::string &reads_params_names_str,
	                                                  const std::string &gtf_filename, const std::string &barcodes_filename)
	{
		std::shared_ptr<Merge::MergeStrategyAbstract> merge_strategy =
				Merge::MergeStrategyFactory::get(this->min_merge_fraction, this->min_genes_before_merge,
												 this->min_genes_after_merge, this->max_merge_edit_distance,
												 merge_tags, barcodes_filename, this->barcode2_length);

		CellsDataContainer container(merge_strategy, Estimator::top_print_size);
		auto bam_processor = BamProcessing::BamProcessorFactory::get(filled_bam, reads_params_names_str,
																	 gtf_filename, this->read_prefix_length);

		auto umig_cells_counts = bam_processor->parse_bam_files(files, bam_output, container);
		container.set_initialized();

		container.merge_and_filter(umig_cells_counts);
		return container;
	}

	Estimator::ss_u_hash_t Estimator::get_umis_per_genes_per_cells_count(const CellsDataContainer &genes_container) const
	{
		ss_u_hash_t result;

		for (auto cell_id : genes_container.filtered_cells())
		{
			std::string cell_barcode = genes_container.cell_barcode(cell_id);
			auto &cell = result[cell_barcode];
			for (auto const &gene : genes_container.cell_genes(cell_id))
			{
				cell[gene.first] = gene.second.size();
			}
		}

		return result;
	}

	Estimator::ss_u_hash_t Estimator::get_reads_per_genes_per_cells_count(const CellsDataContainer &genes_container) const
	{
		ss_u_hash_t result;

		for (auto cell_id : genes_container.filtered_cells())
		{
			std::string cell_barcode = genes_container.cell_barcode(cell_id);
			for (auto const &gene : genes_container.cell_genes(cell_id))
			{
				auto &res_gene = result[cell_barcode][gene.first];
				for (auto const &umi : gene.second)
				{
					res_gene += umi.second;
				}
			}
		}

		return result;
	}
}
