#include "BamController.h"

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"
#include "FilledBamParamsParser.h"
#include "ReadMapParamsParser.h"
#include "FilteringBamProcessor.h"

#include <api/BamReader.h>

namespace Estimation
{
namespace BamProcessing
{
	BamController::BamController(const BamTags &tags, bool filled_bam, const std::string &read_param_filenames,
	                             const std::string &gtf_path)
		: _tags(tags)
		, _filled_bam(filled_bam)
		, _read_param_filenames(read_param_filenames)
		, _gtf_path(gtf_path)
	{}

	void BamController::parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
										CellsDataContainer &container) const
	{
		Tools::trace_time("Start parse bams");

		auto processor = std::shared_ptr<BamProcessorAbstract>(new BamProcessor(container, this->_tags, print_result_bams));
		BamController::process_bam_files(bam_files, print_result_bams, processor);

		Tools::trace_time("Bams parsed");
	}

	void BamController::write_filtered_bam_files(const std::vector<std::string> &bam_files,
	                                             const CellsDataContainer &container) const
	{
		Tools::trace_time("Start write filtered bam");

		auto processor = std::shared_ptr<BamProcessorAbstract>(new FilteringBamProcessor(this->_tags, container));
		BamController::process_bam_files(bam_files, true, processor);

		Tools::trace_time("Filtered bam written");
	}

	void BamController::process_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
										  std::shared_ptr<BamProcessorAbstract> processor) const
	{
		std::shared_ptr<ReadParamsParser> parser = this->get_parser(print_result_bams);

		for (auto const &match_level : processor->container().gene_match_level())
		{
			if (match_level.check(UMI::Mark::HAS_INTRONS) && !parser->has_introns())
			{
				throw std::runtime_error("Genes file should have transcript_id tag or intron records for intron/exon search option");
			}
		}

		for (size_t i = 0; i < bam_files.size(); ++i)
		{
			BamController::parse_bam_file(bam_files[i], processor, parser, bam_files.size() == 1);
			processor->trace_state(bam_files[i]);
		}
	}

	void BamController::parse_bam_file(const std::string &bam_name, std::shared_ptr<BamProcessorAbstract> &processor,
									   std::shared_ptr<ReadParamsParser> &parser, bool trace) const
	{
		using namespace BamTools;

		BamReader reader;
		if (!reader.Open(bam_name))
			throw std::runtime_error("Could not open BAM file: " + bam_name);

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
			catch (std::exception error)
			{
				if (unexpected_chromosome_ids.emplace(alignment.RefID).second)
				{
					L_ERR << "ERROR: can't find chromosome, id: " << alignment.RefID;
				}
				continue;
			}

			processor->inc_reads(); // reads with unknown chromosome are not counted
			if (trace && processor->total_reads_num() % 2000000 == 0)
			{
				processor->trace_state(bam_name);
				break;
			}

			BamController::process_alignment(parser, processor, unexpected_chromosomes, chr_name, alignment);
		}

		reader.Close();
	}

	std::shared_ptr<ReadParamsParser> BamController::get_parser(bool save_read_names) const
	{
		if (this->_filled_bam)
			return std::make_shared<FilledBamParamsParser>(this->_gtf_path, this->_tags);

		if (this->_read_param_filenames != "")
			return std::make_shared<ReadMapParamsParser>(this->_gtf_path, save_read_names, this->_read_param_filenames, this->_tags);

		return std::make_shared<ReadParamsParser>(this->_gtf_path, this->_tags);
	}

	void BamController::process_alignment(std::shared_ptr<ReadParamsParser> parser,
	                                      std::shared_ptr<BamProcessorAbstract> processor,
	                                      std::unordered_set<std::string> &unexpected_chromosomes,
	                                      const std::string &chr_name, const BamTools::BamAlignment &alignment) const
	{
		Tools::ReadParameters read_params;
		if (!parser->get_read_params(alignment, read_params))
		{
			processor->inc_cant_parse_num();
			return;
		}

		std::string gene;
		UMI::Mark mark;
		try
		{
			mark = parser->get_gene(chr_name, alignment, gene);
		}
		catch (Tools::RefGenesContainer::ChrNotFoundException ex)
		{
			if (unexpected_chromosomes.emplace(ex.chr_name).second)
			{
				L_WARN << "WARNING: Can't find chromosome '" << ex.chr_name << "'";
			}
			processor->inc_cant_parse_num();
			return;
		}

		processor->write_alignment(alignment, gene, read_params);
		processor->save_read(read_params.cell_barcode(), chr_name, read_params.umi(), gene, mark);
	}
}
}
