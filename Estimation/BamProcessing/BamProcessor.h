#pragma once
#include <Estimation/CellsDataContainer.h>

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
			const bool print_bam;
			size_t total_intergenic_reads;

		protected:
			virtual std::string get_result_bam_name(const std::string &bam_name) const override;

		public:
			BamProcessor(CellsDataContainer &container, bool print_bam);

			void update_bam(const std::string& bam_file, const BamTools::BamReader &reader) override;

			virtual void trace_state(const std::string& trace_prefix) const override;
			virtual void save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi,
								   const std::string& gene, const CellsDataContainer::Mark &umi_mark) override;
			virtual void write_alignment(BamTools::BamAlignment alignment, const std::string& gene,
								 const Tools::ReadParameters &read_params) override;

			virtual const CellsDataContainer& container() const override;
		};
	}
}