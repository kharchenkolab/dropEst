#pragma once

#include "BamProcessorAbstract.h"

namespace Estimation
{
	class CellsDataContainer;
	namespace BamProcessing
	{
		class FilteringBamProcessor : public BamProcessorAbstract
		{
		public:
			FilteringBamProcessor(const CellsDataContainer &container);

			virtual void trace_state(const std::string& bam_file) const override;
			virtual void save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi, const std::string& gene) override;
			virtual void write_alignment(BamTools::BamAlignment alignment, const std::string& gene,
										 const Tools::ReadParameters &read_params) override;
		};
	}
}
