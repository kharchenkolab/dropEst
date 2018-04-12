#include "ResultsPrinter.h"

#include <RInside.h>

#include <Estimation/CellsDataContainer.h>
#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>
#include <Estimation/Merge/MergeProbabilityValidator.h>

using namespace Rcpp;

namespace Estimation
{
	const size_t ResultsPrinter::top_print_size;

	ResultsPrinter::ResultsPrinter(bool write_matrix, bool reads_output, bool validation_stats, bool umi_correction_info)
		: write_matrix(write_matrix)
		, reads_output(reads_output)
		, validation_stats(validation_stats)
		, umi_correction_info(umi_correction_info)
	{}

	void ResultsPrinter::save_results(const CellsDataContainer &container, const std::string &filename) const
	{
		const std::string list_name = "d";

		if (container.filtered_cells().empty())
		{
			L_WARN << "WARNING: filtered cells are empty. Probably, filtration threshold is too strict"
			       << " or you forgot to run 'merge_and_filter'";
		}
		RInside *R = Tools::init_r();

		this->trace_gene_counts(container);

		L_TRACE << "Compiling diagnostic stats: ";
		auto reads_per_chr_per_cell = this->get_reads_per_chr_per_cell_info(container); // Real cells, all UMIs.
		auto saturation_info = this->get_saturation_analysis_info(container); // Filtered cells, requested UMIs.
		auto mean_reads_per_umi = this->get_mean_reads_per_umi(container); // Real cells, all UMIs.
		auto merge_targets = this->get_merge_targets(container); // All cells.
		IntegerVector aligned_reads_per_cb = wrap(container.get_stat_by_real_cells(Stats::TOTAL_READS_PER_CB)); // Real cells, all UMIs
		IntegerVector aligned_umis_per_cb = wrap(container.get_stat_by_real_cells(Stats::TOTAL_UMIS_PER_CB)); // Real cells, all UMIs
		auto requested_umis_per_cb = this->get_requested_umis_per_cb(container); // Real cells, requested UMIs
		auto requested_reads_per_cb = this->get_requested_umis_per_cb(container, true); // Real cells, requested UMIs
		L_TRACE << "Completed.\n";

		(*R)[list_name] = List::create(
				_["cm"] = this->get_count_matrix(container, true),
				_["cm_raw"] = this->get_count_matrix(container, false),
				_["reads_per_chr_per_cells"] = reads_per_chr_per_cell,
				_["mean_reads_per_umi"] = mean_reads_per_umi,
				_["saturation_info"] = saturation_info,
				_["merge_targets"] = merge_targets,
				_["aligned_reads_per_cell"] = aligned_reads_per_cb,
				_["aligned_umis_per_cell"] = aligned_umis_per_cb,
				_["requested_umis_per_cb"] = requested_umis_per_cb,
				_["requested_reads_per_cb"] = requested_reads_per_cb);

		if (this->umi_correction_info)
		{
			auto const rpupc_list_name = list_name + "_rpupc";
			(*R)[rpupc_list_name] = this->get_reads_per_umi_per_cell(container); // Filtered cells, requested UMIs.
			R->parseEvalQ(list_name + "reads_per_umi_per_cell <- " + rpupc_list_name);
		}

		if (this->validation_stats)
		{
			this->save_validation_stats(list_name, container);
		}

		std::string filename_base = this->extract_filename_base(filename);
		this->save_rds(filename_base, list_name);

		if (this->write_matrix) {
			this->save_mtx(list_name, filename_base);
		}

		L_TRACE << "";
	}

	void ResultsPrinter::save_mtx(const std::string &list_name, const std::string &filename_base) const
	{
		RInside *R = Tools::init_r();
		std::string mtx_file = filename_base + ".mtx";

		L_TRACE << "Writing " + mtx_file + " ...";
		R->parseEvalQ("Matrix::writeMM(" + list_name + "$cm, '" + mtx_file + "')");
		R->parseEvalQ("write.table(colnames(" + list_name + "$cm), '" + filename_base + ".cells.tsv', row.names = F, col.names = F, quote = F)");
		R->parseEvalQ("write.table(rownames(" + list_name + "$cm), '" + filename_base + ".genes.tsv', row.names = F, col.names = F, quote = F)");
		L_TRACE << "Completed.";
	}

	std::string ResultsPrinter::extract_filename_base(const std::string &filename) const
	{
		auto extension_pos = filename.find_last_of('.');
		if (extension_pos != std::string::npos && filename.substr(extension_pos + 1) == "rds")
			return filename.substr(0, extension_pos);

		return filename;
	}

	IntegerMatrix ResultsPrinter::create_matrix(const s_vec_t &col_names, const s_vec_t &row_names,
	                            const i_vec_t &counts)
	{
		IntegerMatrix mat(transpose(IntegerMatrix(int(col_names.size()), int(row_names.size()), counts.begin())));

		colnames(mat) = wrap(col_names);
		rownames(mat) = wrap(row_names);

		return mat;
	}

	List ResultsPrinter::get_saturation_analysis_info(const CellsDataContainer &container) const
	{
		L_TRACE << "Saturation info;";
		i_vec_t reads_by_umig;
		s_vec_t reads_by_umig_cbs;
		s_vec_t reads_by_umig_umis;

		for (size_t cell_id = 0; cell_id < container.total_cells_number(); ++cell_id)
		{
			auto const &cell = container.cell(cell_id);
			if (!cell.is_real())
				continue;

			const auto &cb = cell.barcode();
			for (auto const& gene : cell.requested_reads_per_umi_per_gene(container.gene_match_level()))
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

	DataFrame ResultsPrinter::get_reads_per_chr_per_cell_info(Stats::CellChrStatType stat_type,
	                                                          const CellsDataContainer &container) const
	{
		s_vec_t chr_cell_names, chr_names;
		i_vec_t reads_per_chr_per_cell;
		container.get_stat_by_real_cells(stat_type, chr_cell_names, chr_names, reads_per_chr_per_cell);

		return ResultsPrinter::create_matrix(chr_names, chr_cell_names, reads_per_chr_per_cell);
	}

	List ResultsPrinter::get_reads_per_chr_per_cell_info(const CellsDataContainer &container) const
	{
		L_TRACE << "Reads per chromosome per cell;";
		L_TRACE << "Fill exon results";
		auto exon = this->get_reads_per_chr_per_cell_info(Stats::EXON_READS_PER_CHR_PER_CELL, container);

		L_TRACE << "Fill intron results";
		auto intron = this->get_reads_per_chr_per_cell_info(Stats::INTRON_READS_PER_CHR_PER_CELL, container);

		L_TRACE << "Fill intergenic results";
		auto intergenic = this->get_reads_per_chr_per_cell_info(Stats::INTERGENIC_READS_PER_CHR_PER_CELL, container);

		return List::create(_["Exon"] = exon,
		                    _["Intron"] = intron,
		                    _["Intergenic"] = intergenic);
	}

	SEXP ResultsPrinter::get_count_matrix(const CellsDataContainer &container, bool filtered) const
	{
		if (filtered)
		{
			Tools::trace_time("Compiling filtered count matrix");
		}
		else
		{
			Tools::trace_time("Compiling raw count matrix");
		}

		Tools::init_r();

		s_vec_t gene_names, cell_names;
		SEXP cm;
		if (filtered)
		{
			cm = this->get_count_matrix_filtered(container, container.gene_match_level());
		}
		else
		{
			cm = this->get_count_matrix_raw(container);
		}

		Tools::trace_time("Done");
		return cm;
	}

	void ResultsPrinter::trace_gene_counts(const CellsDataContainer &container) const
	{
		std::unordered_map<std::string, unsigned long> gene_counts_map;
		for (size_t cell_id : container.filtered_cells())
		{
			for (auto const &gm : container.cell(cell_id).genes())
			{
				gene_counts_map[container.gene_indexer().get_value(gm.first)] += gm.second.size();
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
		NumericVector reads_per_umis(container.real_cells_number());
		CharacterVector cell_names(container.real_cells_number());

		size_t out_id = 0;
		for (size_t cell_id = 0; cell_id < container.total_cells_number(); ++cell_id)
		{
			auto const &cell = container.cell(cell_id);
			if (!cell.is_real())
				continue;

			size_t umis_num = 0;
			double reads_per_umi = 0.0;
			for (auto const &gene_rec : cell.genes())
			{
				for (auto const &umi_rec : gene_rec.second.umis())
				{
					reads_per_umi += umi_rec.second.read_count();
				}

				umis_num += gene_rec.second.size();
			}
			reads_per_umi /= umis_num;
			reads_per_umis.at(out_id) = reads_per_umi;
			cell_names.at(out_id) = cell.barcode();
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
		size_t total_genes_num = 0;
		for (auto container_cell_id : container.filtered_cells())
		{
			total_genes_num += container.cell(container_cell_id).requested_genes_num();
		}

		StringIndexer cell_indexer, gene_indexer;
		List reads_per_umi_per_cell(total_genes_num);
		std::vector<unsigned> cell_indexes(total_genes_num), gene_indexes(total_genes_num);

		size_t global_gene_id = 0;
		for (auto container_cell_id : container.filtered_cells())
		{
			auto const &cur_cell = container.cell(container_cell_id);
			unsigned cell_id = cell_indexer.add(cur_cell.barcode());

			for (auto const &gene_rpus : cur_cell.requested_reads_per_umi_per_gene(container.gene_match_level()))
			{
				unsigned gene_id = gene_indexer.add(gene_rpus.first);

				List reads_per_umis(gene_rpus.second.size());
				std::vector<std::string> umis(gene_rpus.second.size());
				size_t umi_ind = 0;

				for (auto const &umi_reads : gene_rpus.second)
				{
					auto mean_quality = cur_cell.at(gene_rpus.first).at(umi_reads.first).mean_quality();
					reads_per_umis[umi_ind] = List::create((unsigned) umi_reads.second, mean_quality);
					umis[umi_ind++] = umi_reads.first;
				}

				reads_per_umis.attr("names") = wrap(umis);

				if (global_gene_id >= total_genes_num)
					throw std::out_of_range("Gene id is to large");

				reads_per_umi_per_cell[global_gene_id] = reads_per_umis;
				cell_indexes[global_gene_id] = cell_id;
				gene_indexes[global_gene_id++] = gene_id;
			}
		}

		return List::create(
				_["cells"] = cell_indexer.values(),
				_["genes"] = gene_indexer.values(),
				_["cell_indexes"] = cell_indexes,
				_["gene_indexes"] = gene_indexes,
				_["reads_per_umi"] = reads_per_umi_per_cell);
	}

	List ResultsPrinter::get_merge_targets(const CellsDataContainer &container) const
	{
		L_TRACE << "Merge targets;";
		std::unordered_map<std::string, std::string> merge_targets;
		for (size_t cell_from_id = 0; cell_from_id < container.total_cells_number(); ++cell_from_id)
		{
			size_t cell_to_id = container.merge_targets()[cell_from_id];
			if (cell_from_id == cell_to_id)
				continue;

			auto const &barcode_from = container.cell(cell_from_id).barcode();
			auto const &barcode_to = container.cell(cell_to_id).barcode();
			merge_targets[barcode_from] = barcode_to;
		}

		return wrap(merge_targets);
	}

	SEXP ResultsPrinter::get_count_matrix_filtered(const CellsDataContainer &container, const UMI::Mark::query_t &query_marks) const
	{
		s_vec_t gene_names, cell_names;
		std::unordered_map<std::string, size_t> gene_ids;
		triplets_vec_t triplets;

		for (size_t column_num = 0; column_num < container.filtered_cells().size(); ++column_num)
		{
			auto const &cur_cell = container.cell(container.filtered_cells()[column_num]);
			cell_names.push_back(cur_cell.barcode());

			for (auto const &umis_per_gene : cur_cell.requested_umis_per_gene(query_marks, this->reads_output))
			{
				auto gene_it = gene_ids.emplace(umis_per_gene.first, gene_ids.size());
				if (gene_it.second)
				{
					gene_names.push_back(umis_per_gene.first);
				}

				size_t row_num = gene_it.first->second;

				triplets.emplace_back(row_num, column_num, umis_per_gene.second);
			}
		}

		L_TRACE << gene_ids.size() << " genes, " << cell_names.size() << " cells.";
		return ResultsPrinter::create_matrix(triplets, gene_ids.size(), cell_names.size(), gene_names, cell_names);
	}

	SEXP ResultsPrinter::get_count_matrix_raw(const CellsDataContainer &container) const
	{
		s_vec_t gene_names, cell_names;
		std::unordered_map<std::string, size_t> gene_ids;
		triplets_vec_t triplets;

		size_t column_num = 0;
		for (size_t cell_id = 0; cell_id < container.total_cells_number(); cell_id++)
		{
			auto const &cell = container.cell(cell_id);
			if (!cell.is_real())
				continue;

			cell_names.push_back(cell.barcode());
			for (auto const &gene : cell.genes())
			{
				auto gene_it = gene_ids.emplace(container.gene_indexer().get_value(gene.first), gene_ids.size());
				if (gene_it.second)
				{
					gene_names.push_back(container.gene_indexer().get_value(gene.first));
				}

				size_t row_num = gene_it.first->second;
				size_t cell_value = gene.second.number_of_umis(this->reads_output);

				triplets.emplace_back(row_num, column_num, cell_value);
			}

			column_num++;
		}

		L_TRACE << gene_ids.size() << " genes, " << cell_names.size() << " cells.";
		return ResultsPrinter::create_matrix(triplets, gene_ids.size(), cell_names.size(), gene_names, cell_names);
	}

	IntegerVector ResultsPrinter::get_requested_umis_per_cb(const CellsDataContainer &container, bool return_reads) const
	{
		IntegerVector requested_umis_per_cb(container.real_cells_number());
		CharacterVector cell_names(container.real_cells_number());

		size_t real_cell_id = 0;
		for (size_t i = 0; i < container.total_cells_number(); ++i)
		{
			auto const &cell = container.cell(i);
			if (!cell.is_real())
				continue;

			cell_names[real_cell_id] = cell.barcode();

			size_t cell_expression = 0;
			if (return_reads)
			{
				for (auto const &gene : cell.requested_umis_per_gene(container.gene_match_level(), true))
				{
					cell_expression += gene.second;
				}
			}
			else
			{
				cell_expression = cell.requested_umis_num();
			}

			requested_umis_per_cb[real_cell_id] = cell_expression;
			real_cell_id++;
		}

		requested_umis_per_cb.attr("names") = cell_names;
		return requested_umis_per_cb;
	}

	SEXP ResultsPrinter::create_matrix(const triplets_vec_t &triplets, size_t total_rows, size_t total_cols,
	                                   const s_vec_t &row_names, const s_vec_t &col_names)
	{
		Eigen::SparseMatrix<unsigned> mat(total_rows, total_cols);
		mat.setFromTriplets(triplets.begin(), triplets.end());

		S4 res(wrap(mat));
		res.slot("Dimnames") = List::create(row_names, col_names);
		return res;
	}

	void ResultsPrinter::save_rds(const std::string &filename_base, const std::string &list_name) const
	{
		RInside *R = Tools::init_r();
		std::string rds_filename = filename_base + ".rds";

		L_TRACE << "";
		Tools::trace_time("Writing R data to " + rds_filename + " ...");
		R->parseEvalQ("saveRDS(" + list_name + ", '" + rds_filename + "')");
		Tools::trace_time("Completed");
	}

	void ResultsPrinter::save_intron_exon_matrices(CellsDataContainer &container, const std::string &filename) const
	{
		const std::string list_name = "matrices";
		RInside *R = Tools::init_r();

		Tools::trace_time("Compiling intron/exon matrices");
		List matrices;
		L_TRACE << "Exon";
		matrices["exon"] = this->get_count_matrix_filtered(container, UMI::Mark::get_by_code("e"));
		L_TRACE << "Intron";
		matrices["intron"] = this->get_count_matrix_filtered(container, UMI::Mark::get_by_code("i"));
		L_TRACE << "Intron/exon spanning";
		matrices["spanning"] = this->get_count_matrix_filtered(container, UMI::Mark::get_by_code("BA"));
		Tools::trace_time("Done");

		(*R)[list_name] = matrices;

		std::string filename_base = this->extract_filename_base(filename);
		this->save_rds(filename_base + ".matrices", list_name);
	}

	List ResultsPrinter::get_merge_validation_info(const std::shared_ptr<Merge::PoissonTargetEstimator> &target_estimator,
	                                               const CellsDataContainer &container, unsigned min_ed, unsigned max_ed,
	                                               size_t cb_pairs_num, unsigned log_period) const
	{
		L_TRACE << "Merge validation;";
		Merge::MergeProbabilityValidator validator(target_estimator);
		validator.run_validation(container, min_ed, max_ed, cb_pairs_num, log_period);

		return List::create(
				_["Probability"]=validator.merge_probs(),
				_["UmisPerCell1"]=validator.umis_per_cell1(),
				_["UmisPerCell2"]=validator.umis_per_cell2(),
				_["EditDistance"]=validator.edit_distances(),
				_["IntersectionSize"]=validator.intersection_size(),
				_["ExpectedIntersectionSize"]=validator.expected_intersection_size()
		);
	}

	void ResultsPrinter::save_validation_stats(const std::string &list_name, const CellsDataContainer &container) const
	{
		const std::string val_info_list_name = list_name + "_val";

		RInside *R = Tools::init_r();

		auto estimator = std::make_shared<Merge::PoissonTargetEstimator>(1, 1);
		estimator->init(container.umi_distribution());

		auto distant = this->get_merge_validation_info(estimator, container, 5, 100, 1000000, 100000); // Real cells, doesn't depend on UMIs
		auto adjacent = this->get_merge_validation_info(estimator, container, 1, 1, 100000, 10000); // Real cells, doesn't depend on UMIs

		(*R)[val_info_list_name] = List::create(_["distant"] = distant, _["adjacent"] = adjacent);
		R->parseEvalQ(list_name + "$merge_validation_info <- " + val_info_list_name);
		R->parseEvalQ("rm(" + val_info_list_name + ")");
	}
}
