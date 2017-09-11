#include "TagsFinderBase.h"

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace TagsSearch
{
	TagsFinderBase::TagsFinderBase(std::shared_ptr<FilesProcessor> files_processor, const boost::property_tree::ptree &processing_config)
		: max_reads(processing_config.get<size_t>("reads_per_out_file", std::numeric_limits<size_t>::max()))
		, min_read_len(processing_config.get<unsigned>("min_align_length", 10))
		, poly_a(processing_config.get<std::string>("poly_a_tail", "AAAAAAAA"))
		, _files_processor(files_processor)
		, _trims_counter()
	{
		srand(time(nullptr));
	}

	void TagsFinderBase::run(bool save_reads_names, bool save_stats)
	{
		L_TRACE << "reading reads ";

		long total_reads_read = 1, parsed_reads = 0;
		std::string file_uid = std::to_string(rand()) + char(rand() % 25 + 'A');

		while (true)
		{
			Tools::ReadParameters params;
			FilesProcessor::FastQRecord r2_record;
			if (!this->parse_fastq_record(r2_record, params))
				break;

			if (total_reads_read % 1000000 == 0)
			{
				L_TRACE << "Total " << total_reads_read << " read (" << parsed_reads << " parsed)";
			}

			++total_reads_read;
			if (params.is_empty() || r2_record.sequence.length() < this->min_read_len)
				continue;

			++parsed_reads;

			std::string read_prefix = "@" + file_uid + std::to_string(total_reads_read);
			if (save_reads_names)
			{
				this->_files_processor->write_read_params(read_prefix, params);
			}

			if (save_stats)
			{
				this->_num_reads_per_cb[params.cell_barcode()]++;
			}

			std::string text = params.encoded_id(read_prefix) + "\n" + r2_record.sequence + "\n" +
					r2_record.description + "\n" + r2_record.quality + "\n";

			bool new_file = this->_files_processor->write(text, this->max_reads);
			if (new_file)
			{
				L_TRACE << "|";
			}
		}

		--total_reads_read;

		L_TRACE << this->results_to_string(total_reads_read);
	}

	std::string TagsFinderBase::results_to_string(long total_reads_read) const
	{
		std::stringstream ss;
		ss << " (" << total_reads_read << " reads)\n"
		   << this->get_additional_stat(total_reads_read) << "\n"
		   << this->_trims_counter.print();

		return ss.str();
	}

	void TagsFinderBase::trim(const std::string &barcodes_tail, std::string &sequence, std::string &quality)
	{
		if (sequence.length() != quality.length())
			throw std::runtime_error("Read has different lengths of sequence and quality string: '" +
											 sequence + "', '" + quality + "'");

		len_t trim_pos = sequence.length();
		// attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
		// RC of UMI+second barcode (up to a length r1_rc_length - spacer_finder parameter)
		std::string rcb = this->rc.rc(barcodes_tail);

		len_t rc_pos = sequence.find(rcb);
		if (rc_pos != std::string::npos)
		{
			trim_pos = rc_pos;
			this->_trims_counter.inc(TrimsCounter::RC);
		}
		else
		{
			// attempt 2: find polyA block
			rc_pos = sequence.find(this->poly_a);
			if (rc_pos != std::string::npos)
			{
				trim_pos = rc_pos;
				this->_trims_counter.inc(TrimsCounter::POLY_A);
			}
		}

		// attempt 3: trim trailing As
		bool a_trim = false;
		len_t skip_count = 0;
		while (trim_pos > 0 && (sequence.at(trim_pos - 1) == 'A' || sequence.at(trim_pos - 1) == 'N'))
		{
			trim_pos--;
			skip_count++;
			a_trim = true;
		}
		if (a_trim)
		{
			this->_trims_counter.inc(TrimsCounter::A_TRIM);
		}

		//attempt 4: apply
		if (sequence.length() != trim_pos)
		{
			sequence = sequence.substr(0, trim_pos);
			quality = quality.substr(0, trim_pos);
		}
		else
		{
			this->_trims_counter.inc(TrimsCounter::NO_TRIM);
		}
	}

	const TagsFinderBase::s_counter_t& TagsFinderBase::num_reads_per_cb() const
	{
		return this->_num_reads_per_cb;
	}
}
