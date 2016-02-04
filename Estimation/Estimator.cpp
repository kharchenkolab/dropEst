#include "Estimator.h"

#include "Tools/log_defs.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

//#include <boost/range/adaptor/reversed.hpp>

using namespace std;

const size_t Estimator::top_print_size = 10;

static bool comp_counters(const pair<string, int> &p1, const pair<string, int> &p2)
{
	return p1.second > p2.second;
}

Estimator::Estimator(const boost::property_tree::ptree &config)
{
	this->read_prefix_length = config.get<size_t>("read_prefix_length");
	this->min_merge_fraction = config.get<double>("min_merge_fraction");
	this->max_merge_edit_distance = config.get<int>("max_merge_edit_distance");
	this->min_genes_after_merge = config.get<int>("min_genes_after_merge");
	this->min_genes_before_merge = config.get<int>("min_genes_before_merge");

	if (this->min_genes_after_merge > 0 && this->min_genes_after_merge < this->min_genes_before_merge)
	{
		this->min_genes_before_merge = this->min_genes_after_merge;
	}
}

IndropResult Estimator::get_results(const names_t &files, bool merge_tags)
{
	GenesContainer container(files, merge_tags, this->read_prefix_length, this->min_merge_fraction,
							 this->min_genes_before_merge, this->min_genes_after_merge,
							 this->max_merge_edit_distance, Estimator::top_print_size);

	ids_t filtered_cells = container.filtered_cells();

	L_TRACE << this->get_cb_top_verbose(container, filtered_cells);

	L_TRACE << filtered_cells.size() << " valid (with >=" << this->min_genes_after_merge << " genes) cells with ";

	s_counter_t gene_counts(this->count_genes(container, filtered_cells));

	L_TRACE << "compiling count matrix ... ";
	names_t cell_names(this->get_unmerged_names(container, filtered_cells));

	names_t gene_names;
	ints_t umis_table;
	this->fill_table(filtered_cells, gene_counts, container, gene_names, umis_table);

	L_TRACE << "Done";

	CountMatrix cm(cell_names, gene_names, umis_table);
	return this->get_indrop_results(cm, container, filtered_cells);
}

Estimator::s_counter_t Estimator::count_genes(const GenesContainer &genes_container, const ids_t &unmerged_cells) const
{
	SIHM counter;
//	for (size_t cell_index: unmerged_cells)
	for (ids_t::const_iterator cell_id_it = unmerged_cells.begin(); cell_id_it != unmerged_cells.end(); ++cell_id_it)
	{
		const GenesContainer::genes_t &cell_genes = genes_container.cell_genes(*cell_id_it);
		for (GenesContainer::genes_t::const_iterator gm_it = cell_genes.begin(); gm_it != cell_genes.end(); ++gm_it)
		{
			counter[gm_it->first] += gm_it->second.size();
		}
	}

	L_TRACE << counter.size() << " genes";

	s_counter_t gene_counts(counter.begin(), counter.end());
	sort(gene_counts.begin(), gene_counts.end(), comp_counters);

	L_TRACE << this->get_genes_top_verbose(gene_counts);
	return gene_counts;
}

Estimator::names_t Estimator::get_unmerged_names(const GenesContainer &genes_container, const ids_t &unmerged_cells) const
{
	names_t cell_names;
	cell_names.reserve(unmerged_cells.size());
	for (ids_t::const_iterator cell_id_it = unmerged_cells.begin(); cell_id_it != unmerged_cells.end(); ++cell_id_it)
	{
		cell_names.push_back(genes_container.gene_name(*cell_id_it));
	}
	return cell_names;
}

void Estimator::fill_table(const ids_t &unmerged_cells, const s_counter_t &gene_counts,
						   const GenesContainer &genes_container, names_t &gene_names_header, ints_t &umis_table) const
{
	size_t size = unmerged_cells.size() * gene_counts.size();
	umis_table.resize(size);
	gene_names_header.reserve(gene_counts.size());

	for (int i = 0; i < gene_counts.size(); i++)
	{
		const string gene_name = gene_counts[i].first;
		gene_names_header.push_back(gene_name);
		for (size_t j = 0; j < unmerged_cells.size(); j++)
		{
			const GenesContainer::genes_t &cur_genes = genes_container.cell_genes(unmerged_cells[j]);
			GenesContainer::genes_t::const_iterator res = cur_genes.find(gene_name);

			umis_table[(i * unmerged_cells.size()) + j] = res == cur_genes.end() ? 0 : res->second.size();
		}
	}
}

IndropResult Estimator::get_indrop_results(const CountMatrix cm, const GenesContainer &genes_container,
										   const ids_t &unmerged_cells) const
{
	L_TRACE << "compiling diagnostic stats: ";

	doubles_t reads_per_umis(this->get_reads_per_umis(genes_container, unmerged_cells));
	L_TRACE << "reads/UMI";

	ints_t umig_coverage(this->get_umig_coverage(genes_container));
	L_TRACE << "UMIg coverage";

	return IndropResult(cm, genes_container.stats(), reads_per_umis, umig_coverage);
}

Estimator::doubles_t Estimator::get_reads_per_umis(const GenesContainer &genes_container, const ids_t &unmerged_cells) const
{
	doubles_t reads_per_umis;
	for (size_t j = 0; j < unmerged_cells.size(); j++)
	{
		size_t umis_count = 0;
		double reads_per_umi = 0.0;
		const GenesContainer::genes_t &cell_genes = genes_container.cell_genes(unmerged_cells[j]);
		for (GenesContainer::genes_t::const_iterator gene_rec_it = cell_genes.begin();
			 gene_rec_it != cell_genes.end(); ++gene_rec_it)
		{
			for (GenesContainer::s_i_hash_t::const_iterator umi_rec_it = gene_rec_it->second.begin();
				 umi_rec_it != gene_rec_it->second.end(); ++umi_rec_it)
			{
				reads_per_umi += umi_rec_it->second;
				umis_count++;
			}
		}
		reads_per_umi /= umis_count;
		reads_per_umis.push_back(reads_per_umi);
	}

	return reads_per_umis;
}

Estimator::ints_t Estimator::get_umig_coverage(const GenesContainer &genes_container) const
{
	ints_t umig_coverage;
	s_set umigs_seen;
//	for (auto const &gene_count : boost::adaptors::reverse(genes_container.cells_genes_counts_sorted()))
	const GenesContainer::i_counter_t &genes_counts = genes_container.cells_genes_counts_sorted();
	for (GenesContainer::i_counter_t::const_reverse_iterator gene_count_it = genes_counts.rbegin();
			gene_count_it != genes_counts.rend(); ++gene_count_it)
	{
		int new_umigs = 0;
		const GenesContainer::genes_t &cell_genes = genes_container.cell_genes(gene_count_it->index);
		for (GenesContainer::genes_t::const_iterator gene_rec_it = cell_genes.begin();
			 gene_rec_it != cell_genes.end(); ++gene_rec_it)
		{
//			for (auto const &umi_rec: gene_rec_it->second)
			for (GenesContainer::s_i_hash_t::const_iterator umi_rec_it = gene_rec_it->second.begin();
					umi_rec_it != gene_rec_it->second.end(); ++umi_rec_it)
			{
				string umig = umi_rec_it->first + gene_rec_it->first;
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


string Estimator::get_cb_top_verbose(const GenesContainer &genes_container, const ids_t &unmerged_cells) const
{
	stringstream ss;
	if (unmerged_cells.size() > 0)
	{
		ss << "top CBs:\n";
		for (size_t i = 0; i < min(unmerged_cells.size(), Estimator::top_print_size); i++)
		{
			ss << genes_container.cell_genes(unmerged_cells[i]).size() << "\t" << genes_container.gene_name(unmerged_cells[i]) << "\n";
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