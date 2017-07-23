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

		this->trace_gene_counts(container);

		L_TRACE << "compiling diagnostic stats: ";
		auto real_cell_names = this->get_real_cell_names(container);
		auto reads_per_chr_per_cell = this->get_reads_per_chr_per_cell_info(container, real_cell_names); // Real cells, all UMIs.
		auto saturation_info = this->get_saturation_analysis_info(container); // Filtered cells, query UMIs.
		auto mean_reads_per_umi = this->get_mean_reads_per_umi(container); // Real cells, all UMIs.
		auto reads_per_umi_per_cell = this->get_reads_per_umi_per_cell(container); // Filtered cells, query UMIs.
		auto merge_targets = this->get_merge_targets(container); // All cells.
		L_TRACE << "Completed.";

		auto count_matrix_filt = this->get_count_matrix(container, true);
		auto count_matrix_raw = this->get_count_matrix(container, false);
		(*R)["d"] = List::create(
				_["cm"] = count_matrix_filt,
				_["cm_raw"] = count_matrix_raw,
				_["reads_per_chr_per_cells"] = reads_per_chr_per_cell,
				_["reads_per_umi_per_cell"] = reads_per_umi_per_cell,
				_["mean_reads_per_umi"] = mean_reads_per_umi,
				_["saturation_info"] = saturation_info,
				_["merge_targets"] = merge_targets,
				_["aligned_reads_per_cell"] = wrap(container.stats().get_raw(Stats::TOTAL_READS_PER_CB)), // All cells. -
				_["aligned_umis_per_cell"] = wrap(container.stats().get_raw(Stats::TOTAL_UMIS_PER_CB)), // All cells. -
				_["fname"] = wrap(filename));

		R->parseEvalQ("dimnames(d$cm$cm) <- list(d$cm$gene_names, d$cm$cell_names); d$cm <- d$cm$cm");
		R->parseEvalQ("dimnames(d$cm_raw$cm) <- list(d$cm_raw$gene_names, d$cm_raw$cell_names); d$cm_raw <- d$cm_raw$cm");
		Tools::trace_time("Writing R data to " + filename + " ...");

		R->parseEvalQ("saveRDS(d, '" + filename + "')");
		Tools::trace_time("Completed");

		if (this->write_matrix) {
			std::string file_trimmed = filename.substr(0, filename.find_last_of("."));
			std::string mtx_file(file_trimmed + ".mtx");
			L_TRACE << "Writing " + mtx_file + " ...";
			R->parseEvalQ("Matrix::writeMM(d$cm, '" + mtx_file + "')");
			R->parseEvalQ("write.table(colnames(d$cm), '" + file_trimmed + ".cells.tsv', row.names = F, col.names = F, quote = F)");
			R->parseEvalQ("write.table(rownames(d$cm), '" + file_trimmed + ".genes.tsv', row.names = F, col.names = F, quote = F)");
			L_TRACE << "Completed.";
		}
	}

	IntegerMatrix ResultsPrinter::create_matrix(const s_vec_t &col_names, const s_vec_t &row_names,
	                            const l_vec_t &counts)
	{
		IntegerMatrix mat(transpose(IntegerMatrix(int(col_names.size()), int(row_names.size()), counts.begin())));

		colnames(mat) = wrap(col_names);
		rownames(mat) = wrap(row_names);

		return mat;
	}

	ResultsPrinter::ResultsPrinter(bool write_matrix, bool reads_output)
		: write_matrix(write_matrix)
		, reads_output(reads_output)
	{}

	List ResultsPrinter::get_saturation_analysis_info(const CellsDataContainer &container) const
	{
		L_TRACE << "Saturation info;";
		l_vec_t reads_by_umig;
		s_vec_t reads_by_umig_cbs;
		s_vec_t reads_by_umig_umis;

		for (size_t cell_id = 0; cell_id < container.cell_barcodes_raw().size(); ++cell_id)
		{
			if (!container.is_cell_real(cell_id))
				continue;

			const auto &cb = container.cell_barcode(cell_id);
			for (auto const& gene : container.query_read_per_umi_per_gene(cell_id))
			{
				for (auto const &reads_per_umi : gene.second)
				{
					reads_by_umig_cbs.push_back(cb);
					reads_by_umig_umis.push_back(reads_per_umi.first);
					reads_by_umig.push_back(reads_per_umi.second);
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
		container.stats().get_filtered(stat_type, cell_names, chr_cell_names, chr_names, reads_per_chr_per_cell);

		return ResultsPrinter::create_matrix(chr_names, chr_cell_names, reads_per_chr_per_cell);
	}

	List ResultsPrinter::get_reads_per_chr_per_cell_info(const CellsDataContainer &container,
	                                                     const ResultsPrinter::s_vec_t &cell_names) const
	{
		L_TRACE << "Reads per chromosome per cell;";
		L_TRACE << "Fill exon results";
		auto exon = this->get_reads_per_chr_per_cell_info(Stats::EXON_READS_PER_CHR_PER_CELL, container, cell_names);

		L_TRACE << "Fill intron results";
		auto intron = this->get_reads_per_chr_per_cell_info(Stats::INTRON_READS_PER_CHR_PER_CELL, container, cell_names);

		L_TRACE << "Fill intergenic results";
		auto intergenic = this->get_reads_per_chr_per_cell_info(Stats::INTERGENIC_READS_PER_CHR_PER_CELL, container,
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

	ResultsPrinter::s_vec_t ResultsPrinter::get_real_cell_names(const CellsDataContainer &container) const
	{
		s_vec_t cell_names;
		cell_names.reserve(container.number_of_real_cells());
		for (size_t cell_id = 0; cell_id < container.cell_barcodes_raw().size(); ++cell_id)
		{
			if (!container.is_cell_real(cell_id))
				continue;

			cell_names.push_back(container.cell_barcode(cell_id));
		}

		return cell_names;
	}

	List ResultsPrinter::get_count_matrix(const CellsDataContainer &container, bool filtered) const
	{
		if (filtered)
		{
			Tools::trace_time("Compiling filtered count matrix");
		}
		else
		{
			Tools::trace_time("Compiling raw count matrix");
		}

		RInside *R = Tools::init_r();
		R->parseEvalQ("library(Matrix)");

		s_vec_t gene_names, cell_names;
		arma::sp_umat cm;
		if (filtered)
		{
			cm = this->get_count_matrix_filtered(container, gene_names, cell_names);
		}
		else
		{
			cm = this->get_count_matrix_raw(container, gene_names, cell_names);
		}


		Tools::trace_time("Done");
		return List::create(_["cm"] = cm, _["gene_names"] = wrap(gene_names), _["cell_names"] = wrap(cell_names));
	}

	void ResultsPrinter::trace_gene_counts(const CellsDataContainer &genes_container) const
	{
		std::unordered_map<std::string, unsigned long> gene_counts_map;
		for (size_t cell_id : genes_container.filtered_cells())
		{
			for (auto const &gm : genes_container.cell_genes(cell_id))
			{
				gene_counts_map[gm.first] += gm.second.size();
			}
		}

		L_TRACE << "\n" << gene_counts_map.size() << " genes";

		using sip_t = std::pair<std::string, unsigned long>;
		std::vector<sip_t> gene_counts(gene_counts_map.begin(), gene_counts_map.end());
		std::sort(gene_counts.begin(), gene_counts.end(),
		          [](const sip_t &p1, const sip_t &p2) { return p1.second > p2.second; });

		std::ostringstream ss;
		ss << "top genes:\n";
		for (size_t i = 0; i < std::min(gene_counts.size(), ResultsPrinter::top_print_size); i++)
		{
			ss << gene_counts[i].first << '\t' << gene_counts[i].second << "\n";
		}

		L_TRACE << ss.str();
	}

	NumericVector ResultsPrinter::get_mean_reads_per_umi(const CellsDataContainer &container) const
	{
		L_TRACE << "Mean reads per UMI;";
		NumericVector reads_per_umis(container.number_of_real_cells());
		CharacterVector cell_names(container.number_of_real_cells());

		size_t out_id = 0;
		for (size_t cell_id = 0; cell_id < container.cell_barcodes_raw().size(); ++cell_id)
		{
			if (!container.is_cell_real(cell_id))
				continue;

			size_t umis_num = 0;
			double reads_per_umi = 0.0;
			for (auto const &gene_rec : container.cell_genes(cell_id))
			{
				for (auto const &umi_rec : gene_rec.second)
				{
					reads_per_umi += umi_rec.second.read_count;
				}

				umis_num += gene_rec.second.size();
			}
			reads_per_umi /= umis_num;
			reads_per_umis.at(out_id) = reads_per_umi;
			cell_names.at(out_id) = container.cell_barcode(cell_id);
			out_id++;
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
			for (auto const &gene_rpus : container.query_read_per_umi_per_gene(cell_id))
			{
				auto &out_gene_umis = cell_reads_p_umigs[gene_rpus.first];
				for (auto const &umi_reads : gene_rpus.second)
				{
					out_gene_umis[umi_reads.first] = (unsigned) umi_reads.second;
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

	arma::sp_umat ResultsPrinter::get_count_matrix_filtered(const CellsDataContainer &container,
	                                                        s_vec_t &gene_names, s_vec_t &cell_names) const
	{
		std::unordered_map<std::string, size_t> gene_ids;
		for (size_t i = 0; i < container.filtered_cells().size(); i++)
		{
			size_t cur_cell_id = container.filtered_cells()[i];
			cell_names.push_back(container.cell_barcode(cur_cell_id));
			for (auto const &umis_per_gene : container.query_umis_per_gene(cur_cell_id, false))
			{
				if (gene_ids.emplace(umis_per_gene.first, gene_ids.size()).second)
				{
					gene_names.push_back(umis_per_gene.first);
				}
			}
		}

		arma::sp_umat cm(arma::uword(gene_ids.size()), arma::uword(container.filtered_cells().size()));
		for (size_t col = 0; col < container.filtered_cells().size(); col++)
		{
			for (auto const &umis_per_gene : container.query_umis_per_gene(container.filtered_cells()[col], this->reads_output))
			{
				size_t row = gene_ids.at(umis_per_gene.first);
				cm(row, col) = arma::uword(umis_per_gene.second);
			}
		}
		return cm;
	}

	arma::sp_umat ResultsPrinter::get_count_matrix_raw(const CellsDataContainer &container,
	                                                   s_vec_t &gene_names, s_vec_t &cell_names) const
	{
		std::unordered_map<std::string, size_t> gene_ids;
		for (size_t cell_id = 0; cell_id < container.cell_barcodes_raw().size(); cell_id++)
		{
			if (!container.is_cell_real(cell_id))
				continue;

			cell_names.push_back(container.cell_barcode(cell_id));
			for (auto const &gene : container.cell_genes(cell_id))
			{
				if (gene_ids.emplace(gene.first, gene_ids.size()).second)
				{
					gene_names.push_back(gene.first);
				}
			}
		}

		size_t column_num = 0;
		arma::sp_umat cm(arma::uword(gene_ids.size()), arma::uword(container.number_of_real_cells()));
		for (size_t cell_id = 0; cell_id < container.cell_barcodes_raw().size(); cell_id++)
		{
			if (!container.is_cell_real(cell_id))
				continue;

			for (auto const &gene : container.cell_genes(cell_id))
			{
				size_t row = gene_ids.at(gene.first);
				size_t cell_value = 0;

				if (this->reads_output)
				{
					for (auto const &umi : gene.second)
					{
						cell_value += umi.second.read_count;
					}
				}
				else
				{
					cell_value = gene.second.size();
				}

				cm(row, column_num) = arma::uword(cell_value);
			}

			column_num++;
		}

		return cm;
	}
}
