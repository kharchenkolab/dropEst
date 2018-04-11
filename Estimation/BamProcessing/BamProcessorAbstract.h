#pragma once

#include <Estimation/CellsDataContainer.h>
#include <Estimation/BamProcessing/BamTags.h>
#include <Estimation/UMI.h>
#include <Estimation/ReadInfo.h>

#include <cstdlib>
#include <string>
#include <api/BamWriter.h>
#include <api/BamReader.h>

namespace Tools
{
	class ReadParameters;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class BamProcessorAbstract
		{
		private:
			size_t _total_reads_num;
			size_t _cant_parse_reads_num;
			size_t _low_quality_reads_num;

			const BamTags _tags;

			BamTools::BamWriter _writer;

		protected:
			void save_alignment(BamTools::BamAlignment alignment, const ReadInfo &read_info_raw,
			                    const Tools::ReadParameters &corrected_params);
			virtual std::string get_result_bam_name(const std::string &bam_name) const = 0;

		public:
			explicit BamProcessorAbstract(const BamTags &tags_info);
			virtual ~BamProcessorAbstract();

			size_t total_reads_num() const;
			size_t cant_parse_reads_num() const;
			size_t low_quality_reads_num() const;
			void inc_reads();
			void inc_cant_parse_num();
			void inc_low_quality_num();
			virtual void update_bam(const std::string& bam_file, const BamTools::BamReader &reader);

			virtual void trace_state(const std::string& trace_prefix) const = 0;
			virtual void save_read(const ReadInfo &read_info) = 0;
			virtual void write_alignment(BamTools::BamAlignment alignment, const ReadInfo &read_info) = 0;

			virtual const CellsDataContainer& container() const = 0;
		};
	}
}