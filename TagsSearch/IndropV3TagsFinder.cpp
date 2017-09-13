#include "IndropV3TagsFinder.h"

namespace TagsSearch
{
	IndropV3TagsFinder::IndropV3TagsFinder(const std::string &barcode1_fastq_name, const std::string &barcode2_fastq_name,
	                                       const std::string &gene_fastq_name,
	                                       const boost::property_tree::ptree &barcodes_config,
	                                       const boost::property_tree::ptree &processing_config,
	                                       TextWriter &&writer, bool save_stats)
		: TagsFinderBase(processing_config, std::move(writer), save_stats)
		, barcode1_length(barcodes_config.get<size_t>("barcode1_length"))
		, barcode2_length(barcodes_config.get<size_t>("barcode2_length"))
		, umi_length(barcodes_config.get<size_t>("umi_length"))
		, trim_tail_length(std::min(barcodes_config.get<size_t>("r1_rc_length"), barcode2_length + umi_length))
		, _barcode1_reader(barcode1_fastq_name)
		, _barcode2_reader(barcode2_fastq_name)
		, _gene_reader(gene_fastq_name)
	{}

	bool IndropV3TagsFinder::parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		auto cb1_rec = this->_barcode1_reader.get_next_record();
		if (cb1_rec.id.empty())
			return false;

		auto cb2_rec = this->_barcode2_reader.get_next_record();
		if (cb2_rec.id.empty())
			throw std::runtime_error("File '" + this->_barcode2_reader.filename() + "', read '" + cb2_rec.id + "': fastq ended prematurely!");

		record = this->_gene_reader.get_next_record();
		if (record.id.empty())
			throw std::runtime_error("File '" + this->_gene_reader.filename() + "', read '" + record.id + "': fastq ended prematurely!");

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
		return this->_counter.print(total_reads_read);
	}

	void IndropV3TagsFinder::skip_records_row()
	{
		this->_barcode1_reader.get_next_record();
		this->_barcode2_reader.get_next_record();
		this->_gene_reader.get_next_record();
	}
}
