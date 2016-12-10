#include <Tools/Logs.h>
#include "SpacerTagsFinder.h"

namespace TagsSearch
{
	SpacerTagsFinder::SpacerTagsFinder(std::shared_ptr<FilesProcessor> files_processor,
									   const boost::property_tree::ptree &barcodes_config,
									   const boost::property_tree::ptree &config)
			: TagsFinderBase(files_processor, config)
			, spacer_finder(barcodes_config)
	{}

	Tools::ReadParameters SpacerTagsFinder::parse_and_trim(const std::string &r1_seq, const std::string &r2_id,
														   std::string &r2_seq, std::string &r2_quality_str)
	{
		if (r2_seq.length() != r2_quality_str.length())
			throw std::runtime_error("Read " + r2_id + " have different length of sequence and quality string: '" +
									 r2_seq + "', '" + r2_quality_str + "'");

		L_DEBUG << r1_seq << ":";

		auto spacer_pos = this->spacer_finder.find_spacer(r1_seq);
		if (spacer_pos.first == SpacerFinder::ERR_CODE)
			return Tools::ReadParameters();

		std::string cell_barcode = this->spacer_finder.parse_cell_barcode(r1_seq, spacer_pos.first, spacer_pos.second);
		L_DEBUG << "-- cell barcode: " << cell_barcode << " (" << cell_barcode.length() << "nt)";

		std::string umi_barcode = this->spacer_finder.parse_umi_barcode(r1_seq, spacer_pos.second);
		if (umi_barcode.length() == 0)
		{
			umi_barcode = "N";
		}

		L_DEBUG << "-- umi barcode: " << umi_barcode;
		L_DEBUG << "R2: " << r2_seq;

		std::string barcodes_tail = this->spacer_finder.parse_r1_rc(r1_seq, spacer_pos.second);
		this->trim(barcodes_tail, r2_seq, r2_quality_str);

		L_DEBUG << " trimmed:" << r2_seq;

		return Tools::ReadParameters(r2_id, cell_barcode, umi_barcode);
	}

	bool SpacerTagsFinder::parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		auto f1_rec = this->files_processor->get_fastq_record(0);
		if (f1_rec.id.empty())
			return false;

		record = this->files_processor->get_fastq_record(1);
		if (record.id.empty())
			throw std::runtime_error("File '" + this->files_processor->filename(1) + "', read '" + record.id + "': fastq ended prematurely!");

		read_params = this->parse_and_trim(f1_rec.sequence, record.id, record.sequence, record.quality);
		return true;
	}

	std::string SpacerTagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->spacer_finder.get_outcomes_counter().print(total_reads_read);
	}
}