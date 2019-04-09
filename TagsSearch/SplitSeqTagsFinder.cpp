#include "SplitSeqTagsFinder.h"

#include <Tools/Logs.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

namespace TagsSearch
{
	SplitSeqTagsFinder::SplitSeqTagsFinder(const std::vector<std::string> &fastq_filenames,
	                                       const boost::property_tree::ptree &barcode_config,
	                                       const boost::property_tree::ptree &processing_config,
	                                       const std::shared_ptr<ConcurrentGzWriter> &writer,
	                                       bool save_stats, bool save_read_params)
		: TagsFinderBase(fastq_filenames, processing_config, writer, save_stats, save_read_params)
		, _barcode_starts(Tools::parse_vec_from_string(barcode_config.get<std::string>("barcode_starts")))
		, _barcode_lengths(Tools::parse_vec_from_string(barcode_config.get<std::string>("barcode_lengths")))
		, _umi_start(barcode_config.get<size_t>("umi_start"))
		, _umi_length(barcode_config.get<size_t>("umi_length"))
		, _min_barcode_read_length(0)
		, _short_seq_read_num(0)
	{
		if (this->_barcode_starts.size() != this->_barcode_lengths.size())
			throw std::runtime_error("barcode_starts and barcode_lengths must have the same number of values");

		if (this->_barcode_starts.empty())
			throw std::runtime_error("You must specify at least one value for barcode_starts");

		for (size_t i = 0; i < this->_barcode_starts.size(); ++i)
		{
			this->_min_barcode_read_length = std::max(this->_min_barcode_read_length, size_t(this->_barcode_starts.at(i) + this->_barcode_lengths.at(i)));
		}
	}

	bool SplitSeqTagsFinder::parse_fastq_records(FastQReader::FastQRecord &gene_record,
												 Tools::ReadParameters &read_params)
	{
		FastQReader::FastQRecord barcode_record;
		if (!this->fastq_reader(0).get_next_record(barcode_record))
			return false;

		if (!this->fastq_reader(1).get_next_record(gene_record))
			throw std::runtime_error("File '" + this->fastq_reader(1).filename() + "', read '" + barcode_record.id + "': fastq ended prematurely!");

		if (barcode_record.sequence.length() < this->_min_barcode_read_length)
		{
			this->_short_seq_read_num++;
			return false;
		}

		std::string cell_barcode = this->parse_cell_barcode(barcode_record.sequence);
		std::string cell_barcode_quality = this->parse_cell_barcode(barcode_record.quality);

		std::string umi = this->parse_umi_barcode(barcode_record.sequence);
		std::string umi_quality = this->parse_umi_barcode(barcode_record.quality);
		read_params = Tools::ReadParameters(cell_barcode, umi, cell_barcode_quality, umi_quality, this->_barcode_phred_threshold);

		return true;
	}

	std::string SplitSeqTagsFinder::parse_cell_barcode(const std::string &sequence)
	{
		std::string cb;
		for (size_t i = 0; i < this->_barcode_starts.size(); ++i)
		{
			cb += sequence.substr(this->_barcode_starts.at(i), this->_barcode_lengths.at(i));
		}

		return cb;
	}

	std::string SplitSeqTagsFinder::parse_umi_barcode(const std::string &sequence)
	{
		return sequence.substr(this->_umi_start, this->_umi_length);
	}

	std::string SplitSeqTagsFinder::get_additional_stat(long total_reads_read) const
	{
		return "Total " + std::to_string(this->_short_seq_read_num) + " short reads (" +
			std::to_string(100 * double(this->_short_seq_read_num) / total_reads_read) + "%)";
	}
}
