#include <Tools/Logs.h>
#include "SpacerTagsFinder.h"

namespace TagsSearch
{
	SpacerTagsFinder::SpacerTagsFinder(std::shared_ptr<FilesProcessor> files_processor,
									   const boost::property_tree::ptree &spacer_config,
									   const boost::property_tree::ptree &config)
			: TagsFinderBase(files_processor, config)
			, spacer_finder(spacer_config)
	{}

	Tools::ReadParameters SpacerTagsFinder::parse(const std::string &r1_seq, const std::string &r1_quality,
		                                              const SpacerFinder::spacer_pos_t &spacer_pos)
	{
		L_DEBUG << r1_seq << ":";

		std::string cell_barcode = this->spacer_finder.parse_cell_barcode(r1_seq, spacer_pos.first, spacer_pos.second);
		std::string cell_barcode_quality = this->spacer_finder.parse_cell_barcode(r1_quality, spacer_pos.first, spacer_pos.second);
		L_DEBUG << "-- cell barcode: " << cell_barcode << " (" << cell_barcode.length() << "nt)";

		std::string umi_barcode = this->spacer_finder.parse_umi_barcode(r1_seq, spacer_pos.second);
		std::string umi_barcode_quality = this->spacer_finder.parse_umi_barcode(r1_quality, spacer_pos.second);
		if (umi_barcode.length() == 0)
		{
			umi_barcode = "N";
			umi_barcode_quality = "@";
		}

		L_DEBUG << "-- umi barcode: " << umi_barcode;

		return Tools::ReadParameters(cell_barcode, umi_barcode, cell_barcode_quality, umi_barcode_quality);
	}

	bool SpacerTagsFinder::parse_fastq_record(FilesProcessor::FastQRecord &gene_record, Tools::ReadParameters &read_params)
	{
		auto barcodes_record = this->_files_processor->get_fastq_record(0);
		if (barcodes_record.id.empty())
			return false;

		gene_record = this->_files_processor->get_fastq_record(1);
		if (gene_record.id.empty())
			throw std::runtime_error("File '" + this->_files_processor->filename(1) + "', read '" + gene_record.id + "': fastq ended prematurely!");

		auto spacer_pos = this->spacer_finder.find_spacer(barcodes_record.sequence);
		if (spacer_pos.first == SpacerFinder::ERR_CODE)
			return true;

		read_params = this->parse(barcodes_record.sequence, barcodes_record.quality, spacer_pos);

		std::string barcodes_tail = this->spacer_finder.parse_r1_rc(barcodes_record.sequence, spacer_pos.second);
		this->trim(barcodes_tail, gene_record.sequence, gene_record.quality);
		return true;
	}

	std::string SpacerTagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->spacer_finder.get_outcomes_counter().print(total_reads_read);
	}
}
