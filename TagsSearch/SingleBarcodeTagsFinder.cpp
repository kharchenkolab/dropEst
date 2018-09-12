#include "SingleBarcodeTagsFinder.h"

namespace TagsSearch
{
	SingleBarcodeTagsFinder::SingleBarcodeTagsFinder(const std::vector<std::string> &fastq_filenames,
	                                                 const boost::property_tree::ptree &barcodes_config,
	                                                 const boost::property_tree::ptree &processing_config,
	                                                 const std::shared_ptr<ConcurrentGzWriter> &writer, bool save_stats, bool save_read_params)
		: TagsFinderBase(fastq_filenames, processing_config, writer, save_stats, save_read_params)
		, barcode_length(barcodes_config.get<size_t>("barcode_length"))
		, umi_length(barcodes_config.get<size_t>("umi_length"))
		, trim_tail_length(std::min(barcodes_config.get<size_t>("r1_rc_length"), barcode_length + umi_length))
	{}

	bool SingleBarcodeTagsFinder::parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		FastQReader::FastQRecord cb_rec;
		if (!this->fastq_reader(0).get_next_record(cb_rec))
			return false;

		if (!this->fastq_reader(1).get_next_record(record))
			throw std::runtime_error("File '" + this->fastq_reader(2).filename() + "', read '" + cb_rec.id + "': fastq ended prematurely!");

		if (cb_rec.sequence.length() < this->barcode_length)
		{
			read_params = Tools::ReadParameters();
			return true;
		}

		std::string cb = this->parse_cb(cb_rec.sequence);
		std::string umi = this->parse_umi(cb_rec.sequence);
		std::string cb_quality = this->parse_cb(cb_rec.quality);
		std::string umi_quality = this->parse_umi(cb_rec.quality);

		if (this->trim_tail_length != 0)
		{
			std::string tail = cb_rec.sequence.substr(this->barcode_length + this->umi_length - this->trim_tail_length, this->trim_tail_length);
			this->trim_poly_a(tail, record.sequence, record.quality);
		}

		read_params = Tools::ReadParameters(cb, umi, cb_quality, umi_quality, this->_barcode_phred_threshold);
		return true;
	}

	std::string SingleBarcodeTagsFinder::parse_umi(const std::string &cb_seq) const
	{
		return cb_seq.substr(this->barcode_length, this->umi_length);
	}

	std::string SingleBarcodeTagsFinder::parse_cb(const std::string &cb_seq) const
	{
		return cb_seq.substr(0, this->barcode_length);
	}
}
