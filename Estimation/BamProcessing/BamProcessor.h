#pragma once
#include <Estimation/CellsDataContainer.h>
#include <Estimation/UMI.h>

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "BamProcessorAbstract.h"

namespace Estimation
{
	namespace BamProcessing
	{
		class BamProcessor : public BamProcessorAbstract
		{
		private:
			CellsDataContainer &_container;
			const bool _print_bam;
			size_t _total_intergenic_reads;

		protected:
			std::string get_result_bam_name(const std::string &bam_name) const override;

		public:
			BamProcessor(CellsDataContainer &container, const BamTags &tags, bool print_bam);

			void update_bam(const std::string& bam_file, const BamTools::BamReader &reader) override;

			void trace_state(const std::string& trace_prefix) const override;

			void save_read(const ReadInfo &read_info) override;

			void write_alignment(BamTools::BamAlignment alignment, const std::string& gene,
								 const Tools::ReadParameters &read_params) override;

			const CellsDataContainer& container() const override;
		};
	}
}