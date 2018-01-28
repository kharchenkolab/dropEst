#include "BamController.h"

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"
#include "FilledBamParamsParser.h"
#include "FilteringBamProcessor.h"
#include "ReadMapParamsParser.h"
#include "Estimation/ReadInfo.h"

#include <api/BamReader.h>

namespace Estimation
{
namespace BamProcessing
{
	BamController::BamController(const BamTags &tags, bool filled_bam, const std::string &read_param_filenames,
	                             const std::string &gtf_path, bool gene_in_chromosome_name, int min_barcode_quality)
		: _tags(tags)
		, _filled_bam(filled_bam)
		, _gene_in_chromosome_name(gene_in_chromosome_name)
		, _read_param_filenames(read_param_filenames)
		, _gtf_path(gtf_path)
		, _min_barcode_quality(min_barcode_quality)
	{}

	void BamController::parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
										CellsDataContainer &container) const
	{
		L_TRACE << "";
		Tools::trace_time("Start parse bams");

		auto processor = std::shared_ptr<BamProcessorAbstract>(new BamProcessor(container, this->_tags, print_result_bams));
		this->process_bam_files(bam_files, processor);

		Tools::trace_time("Bams parsed");
	}

	void BamController::write_filtered_bam_files(const std::vector<std::string> &bam_files,
	                                             const CellsDataContainer &container) const
	{
		L_TRACE << "";
		Tools::trace_time("Start write filtered bam");

		auto processor = std::shared_ptr<BamProcessorAbstract>(new FilteringBamProcessor(this->_tags, container));
		this->process_bam_files(bam_files, processor);

		Tools::trace_time("Filtered bam written");
	}

	void BamController::process_bam_files(const std::vector<std::string> &bam_files,
	                                      std::shared_ptr<BamProcessorAbstract> processor) const
	{
		std::shared_ptr<ReadParamsParser> parser = this->get_parser();

		for (auto const &match_level : processor->container().gene_match_level())
		{
			if (match_level.check(UMI::Mark::HAS_INTRONS) && !parser->has_introns())
			{
				throw std::runtime_error("Genes file should have transcript_id tag or intron records for intron/exon search option");
			}
		}

		for (size_t i = 0; i < bam_files.size(); ++i)
		{
			this->parse_bam_file(bam_files[i], processor, parser, bam_files.size() == 1);
			processor->trace_state(bam_files[i]);
		}
	}

	void BamController::parse_bam_file(const std::string &bam_name, std::shared_ptr<BamProcessorAbstract> &processor,
									   std::shared_ptr<ReadParamsParser> &parser, bool trace) const
	{
		using namespace BamTools;

		BamReader reader;
		if (!reader.Open(bam_name))
			throw std::runtime_error("Can't open BAM file: " + bam_name);

		processor->update_bam(bam_name, reader);

		BamAlignment alignment;
		std::unordered_set<std::string> unexpected_chromosomes;
		std::unordered_set<int32_t> unexpected_chromosome_ids;

		while (reader.GetNextAlignment(alignment))
		{
			std::string chr_name;
			try
			{
				chr_name = reader.GetReferenceData().at(alignment.RefID).RefName;
			}
			catch (std::exception &error)
			{
				if (unexpected_chromosome_ids.emplace(alignment.RefID).second)
				{
					L_ERR << "ERROR: can't find chromosome, id: " << alignment.RefID;
				}
				continue;
			}

			processor->inc_reads(); // reads with unknown chromosome are not counted
			if (trace && processor->total_reads_num() % 10000000 == 0)
			{
				processor->trace_state(bam_name);
			}

			this->process_alignment(parser, processor, unexpected_chromosomes, chr_name, alignment);
		}

		reader.Close();
	}

	std::shared_ptr<ReadParamsParser> BamController::get_parser() const
	{
		if (this->_filled_bam)
			return std::make_shared<FilledBamParamsParser>(this->_gtf_path, this->_tags, this->_gene_in_chromosome_name,
			                                               this->_min_barcode_quality);

		if (this->_read_param_filenames != "")
			return std::make_shared<ReadMapParamsParser>(this->_gtf_path, this->_read_param_filenames,
			                                             this->_tags, this->_gene_in_chromosome_name,
			                                             this->_min_barcode_quality);

		return std::make_shared<ReadParamsParser>(this->_gtf_path, this->_tags, this->_gene_in_chromosome_name);
	}

	void BamController::process_alignment(std::shared_ptr<ReadParamsParser> parser,
	                                      std::shared_ptr<BamProcessorAbstract> processor,
	                                      std::unordered_set<std::string> &unexpected_chromosomes,
	                                      const std::string &chr_name, const BamTools::BamAlignment &alignment) const
	{
		if (!alignment.IsMapped() || !alignment.IsPrimaryAlignment())
			return;

		Tools::ReadParameters read_params;
		if (!parser->get_read_params(alignment, read_params))
		{
			processor->inc_cant_parse_num();
			return;
		}

		if (!read_params.pass_quality_threshold())
		{
			processor->inc_low_quality_num();
			return;
		}

		std::string gene;
		UMI::Mark mark;
		try
		{
			mark = parser->get_gene(chr_name, alignment, gene);
		}
		catch (Tools::GeneAnnotation::RefGenesContainer::ChrNotFoundException ex)
		{
			if (unexpected_chromosomes.emplace(ex.chr_name).second)
			{
				L_WARN << "WARNING: Can't find chromosome '" << ex.chr_name << "'";
			}
			processor->inc_cant_parse_num();
			return;
		}

		ReadInfo read_info(read_params, gene, chr_name, mark);
		processor->write_alignment(alignment, gene, read_params); // TODO: save quality
		processor->save_read(read_info);
	}
}
}
