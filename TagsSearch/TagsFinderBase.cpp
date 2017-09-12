#include "TagsFinderBase.h"

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"

namespace TagsSearch
{
	TagsFinderBase::TagsFinderBase(const boost::property_tree::ptree &processing_config, bool save_stats)
		: _save_stats(save_stats)
		, _file_uid(TagsFinderBase::get_file_uid(42)) // TODO: return time(nullptr)
		, _total_reads_read(0)
		, _parsed_reads(0)
		, _file_ended(false)
		, _max_reads(processing_config.get<size_t>("reads_per_out_file", std::numeric_limits<size_t>::max()))
		, _min_read_len(processing_config.get<unsigned>("min_align_length", 10))
		, poly_a(processing_config.get<std::string>("poly_a_tail", "AAAAAAAA"))
		, _trims_counter()
	{}

	bool TagsFinderBase::get_next_record(FastQReader::FastQRecord& record)
	{
		if (this->_file_ended)
			return false;

		Tools::ReadParameters params;
		if (!this->parse_fastq_record(record, params))
		{
			this->_file_ended = true;
			return false;
		}

		if (++this->_total_reads_read % 1000000 == 0)
		{
			L_TRACE << "Total " << this->_total_reads_read << " read (" << this->_parsed_reads << " parsed)";
		}

		if (params.is_empty() || record.sequence.length() < this->_min_read_len)
			return false;

		++this->_parsed_reads;

		std::string read_prefix = "@" + this->_file_uid + std::to_string(this->_total_reads_read + 1); // TODO: remove "+1"
		record.id = params.encoded_id(read_prefix); // TODO: add save_read_names parameter // record.id = read_prefix;

		if (this->_save_stats)
		{
			this->_num_reads_per_cb[params.cell_barcode()]++;
		}

		return true;
	}

	std::string TagsFinderBase::results_to_string() const
	{
		std::stringstream ss;
		ss << " (" << this->_total_reads_read << " reads)\n"
		   << this->get_additional_stat(this->_total_reads_read) << "\n"
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

	bool TagsFinderBase::file_ended() const
	{
		return this->_file_ended;
	}

	std::string TagsFinderBase::get_file_uid(long random_seed)
	{
		if (random_seed < 0)
		{
			random_seed = time(nullptr);
		}

		srand(unsigned(random_seed));
		return std::to_string(rand()) + char(rand() % 25 + 'A');
	}
}
