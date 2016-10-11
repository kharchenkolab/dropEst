#pragma once

#include "Estimation/CellsDataContainer.h"
#include "Tools/RefGenesContainer.h"
#include "BamProcessor.h"

#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>

namespace Tools
{
	class ReadParameters;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class ReadsParamsParser;
		class BamController
		{
		public:
			static const std::string GENE_TAG;
			static const std::string CB_TAG;
			static const std::string UMI_TAG;

		private:
			static void parse_bam_file(const std::string &bam_name, std::shared_ptr<BamProcessor> &processor,
									   std::shared_ptr<ReadsParamsParser> &parser);

			static std::shared_ptr<ReadsParamsParser> get_parser(bool filled_bam, bool save_read_names,
																 const std::string &reads_params_names_str,
																 const std::string &gtf_path);

			static void parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
										bool filled_bam, const std::string &reads_params_names_str,
										const std::string &gtf_path, std::shared_ptr<BamProcessor> processor);

		public:
			static void parse_bam_files(const std::vector<std::string> &bam_files, bool print_result_bams,
										bool filled_bam, const std::string &reads_params_names_str,
										const std::string &gtf_path, CellsDataContainer &container);
		};
	}
}
