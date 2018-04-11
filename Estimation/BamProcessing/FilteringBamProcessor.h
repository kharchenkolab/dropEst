#pragma once

#include "BamProcessorAbstract.h"
#include <Estimation/CellsDataContainer.h>

#include <unordered_map>
#include <string>

namespace Estimation
{
	namespace BamProcessing
	{
		class FilteringBamProcessor : public BamProcessorAbstract
		{
		private:
			using cb_map_t = std::unordered_map<std::string, std::string>;

			cb_map_t merge_cbs;
			bool is_bam_open;
			size_t written_reads;
			size_t _wrong_genes;
			size_t _wrong_umis;

			const CellsDataContainer &_container;

		protected:
			std::string get_result_bam_name(const std::string &bam_name) const override;

		public:
			FilteringBamProcessor(const BamTags &tags, const CellsDataContainer &container);

			void update_bam(const std::string& bam_file, const BamTools::BamReader &reader) override;
			void trace_state(const std::string& bam_file) const override;
			void save_read(const ReadInfo &read_info) override;
			void write_alignment(BamTools::BamAlignment alignment, const ReadInfo &read_info) override;

			const CellsDataContainer& container() const override;
		};
	}
}
