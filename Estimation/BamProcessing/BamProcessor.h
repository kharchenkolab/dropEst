#pragma once
#include <Estimation/CellsDataContainer.h>

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

namespace Tools
{
	class ReadParameters;
}

namespace Estimation
{
	namespace BamProcessing
	{
		class BamProcessor
		{
		private:
			CellsDataContainer &container;
			CellsDataContainer::s_i_map_t cells_ids;
			const bool print_bam;
			size_t total_exonic_reads;
			size_t total_reads;

			BamTools::BamWriter writer;

		protected:
			virtual std::string get_result_bam_name(const std::string &bam_name) const;

		public:
			BamProcessor(CellsDataContainer &container, bool print_bam);
			virtual ~BamProcessor();

			void inc_reads();
			void update_bam(const std::string& bam_file, const BamTools::BamReader &reader);

			virtual void trace_state(const std::string& bam_file) const;
			virtual void save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi, const std::string& gene);
			virtual void write_alignment(BamTools::BamAlignment alignment, const std::string& gene,
								 const Tools::ReadParameters &read_params);
		};
	}
}