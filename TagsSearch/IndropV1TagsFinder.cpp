#include "IndropV1TagsFinder.h"

#include <Tools/Logs.h>

namespace TagsSearch
{
	IndropV1TagsFinder::IndropV1TagsFinder(const std::string &barcode_fastq_name, const std::string &gene_fastq_name,
	                                       const boost::property_tree::ptree &spacer_config,
	                                       const boost::property_tree::ptree &config, bool save_stats)
		: TagsFinderBase(config, save_stats)
		, _spacer_finder(spacer_config)
		, _barcode_reader(barcode_fastq_name)
		, _gene_reader(gene_fastq_name)
	{}

	Tools::ReadParameters IndropV1TagsFinder::parse(const std::string &r1_seq, const std::string &r1_quality,
		                                              const SpacerFinder::spacer_pos_t &spacer_pos)
	{
		std::string cell_barcode = this->_spacer_finder.parse_cell_barcode(r1_seq, spacer_pos.first, spacer_pos.second);
		std::string cell_barcode_quality = this->_spacer_finder.parse_cell_barcode(r1_quality, spacer_pos.first, spacer_pos.second);

		std::string umi_barcode = this->_spacer_finder.parse_umi_barcode(r1_seq, spacer_pos.second);
		std::string umi_barcode_quality = this->_spacer_finder.parse_umi_barcode(r1_quality, spacer_pos.second);

		return Tools::ReadParameters(cell_barcode, umi_barcode, cell_barcode_quality, umi_barcode_quality);
	}

	bool IndropV1TagsFinder::parse_fastq_record(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params)
	{
		auto barcodes_record = this->_barcode_reader.get_next_record();
		if (barcodes_record.id.empty())
			return false;

		gene_record = this->_gene_reader.get_next_record();
		if (gene_record.id.empty())
			throw std::runtime_error("File '" + this->_gene_reader.filename() + "', read '" + gene_record.id + "': fastq ended prematurely!");

		auto spacer_pos = this->_spacer_finder.find_spacer(barcodes_record.sequence);
		if (spacer_pos.first == SpacerFinder::ERR_CODE)
			return true;

		read_params = this->parse(barcodes_record.sequence, barcodes_record.quality, spacer_pos);

		std::string barcodes_tail = this->_spacer_finder.parse_r1_rc(barcodes_record.sequence, spacer_pos.second);
		this->trim(barcodes_tail, gene_record.sequence, gene_record.quality);
		return true;
	}

	std::string IndropV1TagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->_spacer_finder.get_outcomes_counter().print(total_reads_read);
	}
}
