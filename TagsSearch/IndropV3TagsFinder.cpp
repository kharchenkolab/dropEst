#include "IndropV3TagsFinder.h"

namespace TagsSearch
{
	IndropV3TagsFinder::IndropV3TagsFinder(const std::shared_ptr<FilesProcessor> &files_processor,
												 const boost::property_tree::ptree &barcodes_config,
												 const boost::property_tree::ptree &processing_config)
			: TagsFinderBase(files_processor, processing_config)
			, barcode1_length(barcodes_config.get<size_t>("barcode1_length"))
			, barcode2_length(barcodes_config.get<size_t>("barcode2_length"))
			, umi_length(barcodes_config.get<size_t>("umi_length"))
			, trim_tail_length(std::min(barcodes_config.get<size_t>("r1_rc_length"), barcode2_length + umi_length))
	{}

	bool IndropV3TagsFinder::parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		auto cb1_rec = this->_files_processor->get_fastq_record(0);
		if (cb1_rec.id.empty())
			return false;

		auto cb2_rec = this->_files_processor->get_fastq_record(1);
		if (cb2_rec.id.empty())
			throw std::runtime_error("File '" + this->_files_processor->filename(1) + "', read '" + cb2_rec.id + "': fastq ended prematurely!");

		record = this->_files_processor->get_fastq_record(2);
		if (record.id.empty())
			throw std::runtime_error("File '" + this->_files_processor->filename(2) + "', read '" + record.id + "': fastq ended prematurely!");

		if (cb1_rec.sequence.length() < this->barcode1_length)
		{
			this->counter.inc(TwoBarcodesCounter::SHORT_READ1);
			read_params = Tools::ReadParameters();
			return true;
		}

		if (cb2_rec.sequence.length() < this->barcode2_length + this->umi_length)
		{
			this->counter.inc(TwoBarcodesCounter::SHORT_READ2);
			read_params = Tools::ReadParameters();
			return true;
		}

		this->counter.inc(TwoBarcodesCounter::OK);
		std::string cb = this->parse_cb(cb1_rec.sequence, cb2_rec.sequence);
		std::string umi = this->parse_umi(cb2_rec.sequence);
		std::string cb_quality = this->parse_cb(cb1_rec.quality, cb2_rec.quality);
		std::string umi_quality = this->parse_umi(cb2_rec.quality);

		if (this->trim_tail_length != 0)
		{
			std::string tail = cb2_rec.sequence.substr(this->barcode2_length + this->umi_length - this->trim_tail_length, this->trim_tail_length);
			this->trim(tail, record.sequence, record.quality);
		}

		read_params = Tools::ReadParameters(cb, umi, cb_quality, umi_quality);
		return true;
	}

	std::string IndropV3TagsFinder::parse_umi(const std::string &cb2_seq) const
	{
		return cb2_seq.substr(this->barcode2_length, this->umi_length);
	}

	std::string IndropV3TagsFinder::parse_cb(const std::string &cb1_seq, const std::string &cb2_seq) const
	{
		std::string cb1 = cb1_seq.substr(0, this->barcode1_length);
		std::string cb2 = cb2_seq.substr(0, this->barcode2_length);
		return cb1 + cb2;
	}

	std::string IndropV3TagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->counter.print(total_reads_read);
	}
}
