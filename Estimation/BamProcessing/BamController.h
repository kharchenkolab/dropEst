#pragma once

#include "Estimation/CellsDataContainer.h"
#include "Tools/RefGenesContainer.h"
#include "BamProcessorAbstract.h"

#include <string>
#include <vector>

namespace Tools
{
	class ReadParameters;
}

namespace TestEstimator
{
	struct testGeneMatchLevelUmiExclusion;
	struct testGeneMatchLevelUmiExclusion2;
	struct testPseudoAlignersGenes;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadParamsParser;
		class BamController
		{
			friend struct TestEstimator::testGeneMatchLevelUmiExclusion;
			friend struct TestEstimator::testGeneMatchLevelUmiExclusion2;
			friend struct TestEstimator::testPseudoAlignersGenes;

		private:
			const BamTags _tags;
			const bool _filled_bam;
			const bool _gene_in_chromosome_name;
			const std::string _read_param_filenames;
			const std::string _gtf_path;

		private:
			void parse_bam_file(const std::string &bam_name, std::shared_ptr<BamProcessorAbstract> &processor,
			                    std::shared_ptr<ReadParamsParser> &parser, bool trace) const;

			std::shared_ptr<ReadParamsParser> get_parser() const;

			void process_bam_files(const std::vector<std::string> &bam_files,
			                       std::shared_ptr<BamProcessorAbstract> processor) const;

			void process_alignment(std::shared_ptr<ReadParamsParser> parser,
			                       std::shared_ptr<BamProcessorAbstract> processor,
			                       std::unordered_set<std::string> &unexpected_chromosomes,
			                       const std::string &chr_name, const BamTools::BamAlignment &alignment) const;

		public:
			void parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
			                     CellsDataContainer &container) const;

			void write_filtered_bam_files(const std::vector<std::string> &bam_files,
			                              const CellsDataContainer &container) const;

			BamController(const BamTags &tags, bool filled_bam, const std::string &read_param_filenames,
			              const std::string &gtf_path, bool gene_in_chromosome_name);
		};
	}
}
