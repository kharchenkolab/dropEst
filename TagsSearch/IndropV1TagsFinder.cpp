#include "IndropV1TagsFinder.h"

#include <Tools/Logs.h>

namespace TagsSearch
{
	IndropV1TagsFinder::IndropV1TagsFinder(const std::vector<std::string> &fastq_filenames,
	                                       const boost::property_tree::ptree &spacer_config,
	                                       const boost::property_tree::ptree &config,
	                                       const std::shared_ptr<ConcurrentGzWriter> &writer,
	                                       bool save_stats, bool save_read_params)
		: TagsFinderBase(fastq_filenames, config, writer, save_stats, save_read_params)
		, _spacer_finder(spacer_config)
	{}

	Tools::ReadParameters IndropV1TagsFinder::parse(const std::string &r1_seq, const std::string &r1_quality,
		                                              const SpacerFinder::spacer_pos_t &spacer_pos)
	{
		std::string cell_barcode = this->_spacer_finder.parse_cell_barcode(r1_seq, spacer_pos.first, spacer_pos.second);
		std::string cell_barcode_quality = this->_spacer_finder.parse_cell_barcode(r1_quality, spacer_pos.first, spacer_pos.second);

		std::string umi_barcode = this->_spacer_finder.parse_umi_barcode(r1_seq, spacer_pos.second);
		std::string umi_barcode_quality = this->_spacer_finder.parse_umi_barcode(r1_quality, spacer_pos.second);

		return Tools::ReadParameters(cell_barcode, umi_barcode, cell_barcode_quality, umi_barcode_quality, this->_barcode_phred_threshold);
	}

	bool IndropV1TagsFinder::parse_fastq_record(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params)
	{
		FastQReader::FastQRecord barcodes_record;
		if (!this->fastq_reader(0).get_next_record(barcodes_record))
			return false;

		if (!this->fastq_reader(1).get_next_record(gene_record))
			throw std::runtime_error("File '" + this->fastq_reader(1).filename() + "', read '" + barcodes_record.id + "': fastq ended prematurely!");

		auto spacer_pos = this->_spacer_finder.find_spacer(barcodes_record.sequence);
		if (spacer_pos.first == SpacerFinder::ERR_CODE)
			return true;

		read_params = this->parse(barcodes_record.sequence, barcodes_record.quality, spacer_pos);

		std::string barcodes_tail = this->_spacer_finder.parse_r1_rc(barcodes_record.sequence, spacer_pos.second);
		this->trim_poly_a(barcodes_tail, gene_record.sequence, gene_record.quality);
		return true;
	}

	std::string IndropV1TagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->_spacer_finder.get_outcomes_counter().print(total_reads_read);
	}
}
