#include <Tools/Logs.h>
#include "FilteringBamProcessor.h"
#include <Estimation/CellsDataContainer.h>

namespace Estimation
{
	namespace BamProcessing
	{

		FilteringBamProcessor::FilteringBamProcessor(const CellsDataContainer &container)
		{
			auto const &merge_targets = container.merge_targets();
			for (size_t base_cell_id = 0; base_cell_id < merge_targets.size(); ++base_cell_id)
			{
				if (container.is_cell_excluded(base_cell_id))
					continue;

				size_t target_id = merge_targets[base_cell_id];
				if (target_id == base_cell_id)
					continue;
			}
		}

		void FilteringBamProcessor::trace_state(const std::string &bam_file) const
		{
			L_TRACE << bam_file << ": " << this->total_reads() << " reads processed";
		}

		void FilteringBamProcessor::save_read(const std::string &cell_barcode, const std::string &chr_name,
											  const std::string &umi, const std::string &gene)
		{}

		void FilteringBamProcessor::write_alignment(BamTools::BamAlignment alignment, const std::string &gene,
													const Tools::ReadParameters &read_params)
		{
		}
	}
}