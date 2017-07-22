#include "ResultsPrinter.h"

#include <RInside.h>

#include <Estimation/CellsDataContainer.h>
#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>

using namespace Rcpp;

namespace Estimation
{
	const size_t ResultsPrinter::top_print_size;

	void ResultsPrinter::save_results(const CellsDataContainer &container, const std::string &filename) const
	{
		if (container.filtered_cells().empty())
		{
			L_WARN << "WARNING: filtered cells is empty. Probably, filtration threshold is too strict"
			       << " or you forgot to run 'merge_and_filter'";
		}
		RInside *R = Tools::init_r();

		auto gene_counts = this->get_gene_counts_sorted(container);

		L_TRACE << "compiling diagnostic stats: ";
		auto cell_names = this->get_filtered_cell_names(container);
		auto reads_per_chr_per_cell = this->get_reads_per_chr_per_cell_info(container, cell_names);
		auto saturation_info = this->get_saturation_analysis_info(container);
		auto mean_reads_per_umi = this->get_mean_reads_per_umi(container);
		auto reads_per_umi_per_cell = this->get_reads_per_umi_per_cell(container);
		auto merge_targets = this->get_merge_targets(container);
		L_TRACE << "Completed.";

		auto count_matrix = this->get_count_matrix(container, gene_counts);
		(*R)["data"] = List::create(
				_["cm"] = count_matrix,
				_["reads_per_chr_per_cells"] = reads_per_chr_per_cell,
				_["reads_per_umi_per_cell"] = reads_per_umi_per_cell,
				_["mean_reads_per_umi"] = mean_reads_per_umi,
				_["saturation_info"] = saturation_info,
				_["merge_targets"] = merge_targets,
				_["aligned_reads_per_cell"] = wrap(container.stats().get_raw(Stats::TOTAL_READS_PER_CB)),
				_["aligned_umis_per_cell"] = wrap(container.stats().get_raw(Stats::TOTAL_UMIS_PER_CB)),
				_["fname"] = wrap(filename));

		Tools::trace_time("Writing R data to " + filename);

		R->parseEvalQ("saveRDS(data, '" + filename + "')");
		Tools::trace_time("Done");
	}

	IntegerMatrix ResultsPrinter::create_matrix(const s_vec_t &col_names, const s_vec_t &row_names,
	                            const l_vec_t &counts)
	{
		IntegerMatrix mat(transpose(IntegerMatrix(int(col_names.size()), int(row_names.size()), counts.begin())));

		colnames(mat) = wrap(col_names);
		rownames(mat) = wrap(row_names);

		return mat;
	}

	ResultsPrinter::ResultsPrinter(bool text_output, bool reads_output)
		: text_output(text_output)
		, reads_output(reads_output)
	{}

	List ResultsPrinter::get_saturation_analysis_info(const CellsDataContainer &container) const
	{
		L_TRACE << "Saturation info;";
		l_vec_t reads_by_umig;
		s_vec_t reads_by_umig_cbs;
		s_vec_t reads_by_umig_umis;

		for (size_t cell_ind = 0; cell_ind < container.cell_barcodes_raw().size(); ++cell_ind)
		{
			if (container.is_cell_merged(cell_ind))
				continue;

			const auto &cb = container.cell_barcode(cell_ind);
			for (auto const& gene : container.cell_genes(cell_ind))
			{
				for (auto const &umi : gene.second)
				{
					reads_by_umig_cbs.push_back(cb);
					reads_by_umig_umis.push_back(umi.first);
					reads_by_umig.push_back(umi.second.read_count);
				}
			}
		}
		return List::create(
				_["reads"] = wrap(reads_by_umig),
				_["cbs"] = wrap(reads_by_umig_cbs),
				_["umis"] = wrap(reads_by_umig_umis)
		);
	}

	DataFrame ResultsPrinter::get_reads_per_chr_per_cell_info(Stats::CellStrStatType stat_type,
	                                                          const CellsDataContainer &container,
	                                                          const s_vec_t &cell_names) const
	{
		s_vec_t chr_cell_names, chr_names;
		l_vec_t reads_per_chr_per_cell;
		container.stats().get(stat_type, chr_cell_names, chr_names, reads_per_chr_per_cell);

		return ResultsPrinter::create_matrix(chr_names, chr_cell_names, reads_per_chr_per_cell);
	}

	List ResultsPrinter::get_reads_per_chr_per_cell_info(const CellsDataContainer &container,
	                                                          const ResultsPrinter::s_vec_t &cell_names) const
	{
		L_TRACE << "Reads per chromosome per cell;";
		L_TRACE << "Fill exon results";
		auto exon = get_reads_per_chr_per_cell_info(Stats::EXON_READS_PER_CHR_PER_CELL, container, cell_names);

		L_TRACE << "Fill intron results";
		auto intron = get_reads_per_chr_per_cell_info(Stats::INTRON_READS_PER_CHR_PER_CELL, container, cell_names);

		L_TRACE << "Fill intergenic results";
		auto intergenic = get_reads_per_chr_per_cell_info(Stats::INTERGENIC_READS_PER_CHR_PER_CELL, container,
		                                                  cell_names);

		return List::create(_["Exon"] = exon,
		                    _["Intron"] = intron,
		                    _["Intergenic"] = intergenic);
	}

	ResultsPrinter::s_vec_t ResultsPrinter::get_filtered_cell_names(const CellsDataContainer &container) const
	{
		s_vec_t cell_names;
		cell_names.reserve(container.filtered_cells().size());
		for (auto const &cell_id : container.filtered_cells())
		{
			cell_names.push_back(container.cell_barcode(cell_id));
		}
		return cell_names;
	}

	IntegerMatrix
	ResultsPrinter::get_count_matrix(const CellsDataContainer &container, const s_counter_t &gene_counts) const
	{
		Tools::trace_time("Compiling count matrix");
		s_vec_t cell_names, gene_names(gene_counts.size());
		l_vec_t counts(container.filtered_cells().size() * gene_names.size(), 0);
		for (int i = 0; i < gene_counts.size(); ++i)
		{
			gene_names[i] = gene_counts[i].first;
		}
		std::unordered_map<std::string, size_t> gene_ids;

		for (size_t i = 0; i < gene_names.size(); ++i)
		{
			gene_ids[gene_names[i]] = i;
		}

		for (size_t col = 0; col < container.filtered_cells().size(); col++)
		{
			size_t cur_cell_id = container.filtered_cells()[col];
			cell_names.push_back(container.cell_barcode(cur_cell_id));
			for (auto const &gene_it : container.cell_genes(cur_cell_id))
			{
				size_t row = gene_ids[gene_it.first];
				int cell_value = 0;

				if (this->reads_output)
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

				counts[(row * container.filtered_cells().size()) + col] = cell_value;
			}
		}

		Tools::trace_time("Done");
		return ResultsPrinter::create_matrix(cell_names, gene_names, counts);
	}

	ResultsPrinter::s_counter_t ResultsPrinter::get_gene_counts_sorted(const CellsDataContainer &genes_container) const
	{
		std::unordered_map<std::string, unsigned long> gene_counts_map;
		for (size_t cell_index : genes_container.filtered_cells())
		{
			for (auto const &gm : genes_container.cell_genes(cell_index))
			{
				gene_counts_map[gm.first] += gm.second.size();
			}
		}

		L_TRACE << "\n" << gene_counts_map.size() << " genes";

		s_counter_t gene_counts(gene_counts_map.begin(), gene_counts_map.end());
		std::sort(gene_counts.begin(), gene_counts.end(),
		          [](const si_pair_t &p1, const si_pair_t &p2) { return p1.second > p2.second; });

		std::ostringstream ss;
		ss << "top genes:\n";
		for (size_t i = 0; i < std::min(gene_counts.size(), ResultsPrinter::top_print_size); i++)
		{
			ss << gene_counts[i].first << '\t' << gene_counts[i].second << "\n";
		}

		L_TRACE << ss.str();

		return gene_counts;
	}

	NumericVector ResultsPrinter::get_mean_reads_per_umi(const CellsDataContainer &container) const
	{
		L_TRACE << "Mean reads per UMI;";
		NumericVector reads_per_umis(container.filtered_cells().size());
		CharacterVector cell_names(container.filtered_cells().size());

		for (size_t i = 0; i < container.filtered_cells().size(); ++i)
		{
			size_t umis_count = 0;
			double reads_per_umi = 0.0;
			for (auto const &gene_rec : container.cell_genes(container.filtered_cells()[i]))
			{
				for (auto const &umi_rec : gene_rec.second)
				{
					reads_per_umi += umi_rec.second.read_count;
				}

				umis_count += gene_rec.second.size();
			}
			reads_per_umi /= umis_count;
			reads_per_umis[i] = reads_per_umi;
			cell_names[i] = container.cell_barcode(i);
		}

		reads_per_umis.attr("names") = cell_names;
		return reads_per_umis;
	}

	List ResultsPrinter::get_reads_per_umi_per_cell(const CellsDataContainer &container) const
	{
		L_TRACE << "Reads per UMI per gene;";
		using std::unordered_map;
		using std::string;
		unordered_map<string, unordered_map<string, unordered_map<string, unsigned>>> reads_per_umi_per_cell;
		for (auto cell_id : container.filtered_cells())
		{
			auto &cell_reads_p_umigs = reads_per_umi_per_cell[container.cell_barcode(cell_id)];
			for (auto const &gene_umis : container.cell_genes(cell_id))
			{
				auto &out_gene_umis = cell_reads_p_umigs[gene_umis.first];
				for (auto const &umi_reads : gene_umis.second)
				{
					out_gene_umis[umi_reads.first] = (unsigned) umi_reads.second.read_count;
				}
			}
		}

		return wrap(reads_per_umi_per_cell);
	}

	List ResultsPrinter::get_merge_targets(const CellsDataContainer &container) const
	{
		L_TRACE << "Merge targets;";
		std::unordered_map<std::string, std::string> merge_targets;
		for (size_t cell_from_id = 0; cell_from_id < container.cell_barcodes_raw().size(); ++cell_from_id)
		{
			auto const &barcode_from = container.cell_barcodes_raw()[cell_from_id];
			auto const &barcode_to = container.cell_barcodes_raw()[container.merge_targets()[cell_from_id]];
			merge_targets[barcode_from] = barcode_to;
		}

		return wrap(merge_targets);
	}
}
