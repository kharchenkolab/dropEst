#include "TwoBarcodesTagsFinder.h"

namespace TagsSearch
{
	TwoBarcodesTagsFinder::TwoBarcodesTagsFinder(const std::shared_ptr<FilesProcessor> &files_processor,
												 const boost::property_tree::ptree &barcodes_config,
												 const boost::property_tree::ptree &config)
			: TagsFinderBase(files_processor, config)
			, barcode1_length(barcodes_config.get<size_t>("barcode1_length"))
			, barcode2_length(barcodes_config.get<size_t>("barcode2_length"))
			, umi_length(barcodes_config.get<size_t>("umi_length"))
			, trim_tail_length(std::min(barcodes_config.get<size_t>("r1_rc_length"), barcode2_length + umi_length))
	{}

	bool TwoBarcodesTagsFinder::parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		auto cb1_rec = this->files_processor->get_fastq_record(0);
		if (cb1_rec.id.empty())
			return false;

		auto cb2_rec = this->files_processor->get_fastq_record(1);
		if (cb2_rec.id.empty())
			throw std::runtime_error("File '" + this->files_processor->filename(1) + "', read '" + cb2_rec.id + "': fastq ended prematurely!");

		record = this->files_processor->get_fastq_record(2);
		if (record.id.empty())
			throw std::runtime_error("File '" + this->files_processor->filename(2) + "', read '" + record.id + "': fastq ended prematurely!");

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
		std::string cb1 = cb1_rec.sequence.substr(0, this->barcode1_length);
		std::string cb2 = cb2_rec.sequence.substr(0, this->barcode2_length);
		std::string umi = cb2_rec.sequence.substr(this->barcode2_length, this->umi_length);

		std::string tail = cb2_rec.sequence.substr(this->barcode2_length + this->umi_length - this->trim_tail_length, this->trim_tail_length);
		this->trim(tail, record.sequence, record.quality);
		read_params = Tools::ReadParameters(record.id, cb1 + cb2, umi);
		return true;
	}

	std::string TwoBarcodesTagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->counter.print(total_reads_read);
	}
}
