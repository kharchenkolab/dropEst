#include "IndropV3TagsFinder.h"

namespace TagsSearch
{
	IndropV3TagsFinder::IndropV3TagsFinder(const std::vector<std::string> &fastq_filenames,
	                                       const boost::property_tree::ptree &barcodes_config,
	                                       const boost::property_tree::ptree &processing_config,
	                                       const std::shared_ptr<ConcurrentGzWriter> &writer, bool save_stats, bool save_read_params)
		: TagsFinderBase(fastq_filenames, processing_config, writer, save_stats, save_read_params)
		, barcode1_length(barcodes_config.get<size_t>("barcode1_length"))
		, barcode2_length(barcodes_config.get<size_t>("barcode2_length"))
		, umi_length(barcodes_config.get<size_t>("umi_length"))
		, trim_tail_length(std::min(barcodes_config.get<size_t>("r1_rc_length"), barcode2_length + umi_length))
	{}

	bool IndropV3TagsFinder::parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		FastQReader::FastQRecord cb1_rec, cb2_rec;
		if (!this->fastq_reader(0).get_next_record(cb1_rec))
			return false;

		if (!this->fastq_reader(1).get_next_record(cb2_rec))
			throw std::runtime_error("File '" + this->fastq_reader(1).filename() + "', read '" + cb1_rec.id + "': fastq ended prematurely!");

		if (!this->fastq_reader(2).get_next_record(record))
			throw std::runtime_error("File '" + this->fastq_reader(2).filename() + "', read '" + cb1_rec.id + "': fastq ended prematurely!");

		if (cb1_rec.sequence.length() < this->barcode1_length)
		{
			this->_counter.inc(TwoBarcodesCounter::SHORT_READ1);
			read_params = Tools::ReadParameters();
			return true;
		}

		if (cb2_rec.sequence.length() < this->barcode2_length + this->umi_length)
		{
			this->_counter.inc(TwoBarcodesCounter::SHORT_READ2);
			read_params = Tools::ReadParameters();
			return true;
		}

		this->_counter.inc(TwoBarcodesCounter::OK);
		std::string cb = this->parse_cb(cb1_rec.sequence, cb2_rec.sequence);
		std::string umi = this->parse_umi(cb2_rec.sequence);
		std::string cb_quality = this->parse_cb(cb1_rec.quality, cb2_rec.quality);
		std::string umi_quality = this->parse_umi(cb2_rec.quality);

		if (this->trim_tail_length != 0)
		{
			std::string tail = cb2_rec.sequence.substr(this->barcode2_length + this->umi_length - this->trim_tail_length, this->trim_tail_length);
			this->trim_poly_a(tail, record.sequence, record.quality);
		}

		read_params = Tools::ReadParameters(cb, umi, cb_quality, umi_quality, this->_barcode_phred_threshold);
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
		return this->_counter.print(total_reads_read);
	}
}
