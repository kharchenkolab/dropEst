#pragma once

#include <Estimation/CellsDataContainer.h>

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
			BamTools::BamWriter _writer;

		protected:
			void save_alignment(BamTools::BamAlignment alignment, const std::string &name, const std::string &gene,
								const std::string &barcode, const std::string &umi);
			virtual std::string get_result_bam_name(const std::string &bam_name) const = 0;

		public:
			BamProcessorAbstract();
			virtual ~BamProcessorAbstract();

			size_t cant_parse_reads_num() const;
			size_t total_reads_num() const;
			void inc_reads();
			void inc_cant_parse_num();
			virtual void update_bam(const std::string& bam_file, const BamTools::BamReader &reader);

			virtual void trace_state(const std::string& trace_prefix) const = 0;
			virtual void save_read(const std::string& cell_barcode, const std::string& chr_name, const std::string& umi,
			                       const std::string& gene, const CellsDataContainer::Mark &umi_mark) = 0;
			virtual void write_alignment(BamTools::BamAlignment alignment, const std::string& gene,
								 const Tools::ReadParameters &read_params) = 0;

			virtual const CellsDataContainer& container() const = 0;
		};
	}
}