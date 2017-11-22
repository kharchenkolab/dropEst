#include "IClipTagsFinder.h"

#include <Tools/Logs.h>

namespace TagsSearch
{
	IClipTagsFinder::IClipTagsFinder(const std::vector<std::string> &fastq_filenames,
	                                 const boost::property_tree::ptree &barcodes_config,
	                                 const boost::property_tree::ptree &config,
	                                 const std::shared_ptr<ConcurrentGzWriter> &writer,
	                                 bool save_stats, bool save_read_params)
		: TagsFinderBase(fastq_filenames, config, writer, save_stats, save_read_params)
		, _barcode_length(barcodes_config.get<size_t>("barcode_length"))
		, _umi_length(barcodes_config.get<size_t>("umi_length"))
		, _cant_parse_num(0)
	{}

	bool IClipTagsFinder::parse_fastq_record(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params)
	{
		if (!this->fastq_reader(0).get_next_record(gene_record))
			return false;

		if (gene_record.sequence.length() <= this->_umi_length + this->_barcode_length + this->_min_read_len)
		{
			this->_cant_parse_num++;
			read_params = Tools::ReadParameters();
			return true;
		}

		std::string barcode = this->parse_cb(gene_record.sequence);
		std::string barcode_quality = this->parse_cb(gene_record.quality);

		std::string umi = this->parse_umi(gene_record.sequence);
		std::string umi_quality = this->parse_umi(gene_record.quality);

		gene_record.sequence = this->trim_barcodes(gene_record.sequence);
		gene_record.quality = this->trim_barcodes(gene_record.quality);

		read_params = Tools::ReadParameters(barcode, umi, barcode_quality, umi_quality);

		return true;
	}

	std::string IClipTagsFinder::get_additional_stat(long total_reads_read) const
	{
		return "Can't parse: " + std::to_string(100 * double(this->_cant_parse_num) / total_reads_read) + "%";
	}

	std::string IClipTagsFinder::trim_barcodes(const std::string &sequence) const
	{
		return sequence.substr(this->_umi_length + this->_barcode_length);
	}

	std::string IClipTagsFinder::parse_cb(const std::string &sequence) const
	{
		return sequence.substr(this->_umi_length, this->_barcode_length);
	}

	std::string IClipTagsFinder::parse_umi(const std::string &sequence) const
	{
		return sequence.substr(0, this->_umi_length);
	}
}
