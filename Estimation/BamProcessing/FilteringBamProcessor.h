#pragma once

#include "BamProcessorAbstract.h"

#include <unordered_map>
#include <string>

namespace Estimation
{
	class CellsDataContainer;
	namespace BamProcessing
	{
		class FilteringBamProcessor : public BamProcessorAbstract
		{
		private:
			typedef std::unordered_map<std::string, std::string> cb_map_t;

			cb_map_t merge_cbs;
			bool is_bam_open;
			size_t written_reads;

		protected:
			virtual std::string get_result_bam_name(const std::string &bam_name) const override;

		public:
			FilteringBamProcessor(const CellsDataContainer &container);

			virtual void update_bam(const std::string& bam_file, const BamTools::BamReader &reader) override;
			virtual void trace_state(const std::string& bam_file) const override;
			virtual void save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi,
								   const std::string& gene) override;
			virtual void write_alignment(BamTools::BamAlignment alignment, const std::string& gene,
										 const Tools::ReadParameters &read_params) override;
		};
	}
}
