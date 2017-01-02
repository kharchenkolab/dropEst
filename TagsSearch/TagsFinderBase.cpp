#include "TagsFinderBase.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"
#include "Tools/UtilFunctions.h"

using namespace std;

namespace TagsSearch
{
	TagsFinderBase::TagsFinderBase(std::shared_ptr<FilesProcessor> files_processor, const boost::property_tree::ptree &config)
		: max_reads(config.get<size_t>("max_reads", numeric_limits<size_t>::max()))
		, min_read_len(config.get<unsigned>("min_align_length"))
		, poly_a(config.get<string>("poly_a_tail"))
		, files_processor(files_processor)
		, trims_counter()
	{
		srand(time(nullptr));
	}

	void TagsFinderBase::run(bool save_reads_names)
	{
		L_TRACE << "reading reads ";

		long total_reads_read = 1, parsed_reads = 0;
		std::string file_uid = std::to_string(rand()) + char(rand() % 25 + 'A');

		Tools::ReadParameters params;
		FilesProcessor::FastQRecord r2_record;
		while (this->parse_fastq_record(r2_record, params))
		{
			if (total_reads_read % 1000000 == 0)
			{
				L_TRACE << "Total " << total_reads_read << " read (" << parsed_reads << " parsed)";
			}

			++total_reads_read;
			if (params.is_empty() || r2_record.sequence.length() < this->min_read_len)
				continue;

			++parsed_reads;

			string text;
			std::string new_id = "@" + file_uid + std::to_string(total_reads_read);
			if (save_reads_names)
			{
				text = new_id;
				this->files_processor->write_read_params(new_id, params);
			}
			else
			{
				text = params.to_monolithic_string(new_id);
			}

			text += "\n" + r2_record.sequence + "\n" + r2_record.description + "\n" + r2_record.quality + "\n";

			bool new_file = this->files_processor->write(text, this->max_reads);
			if (new_file)
			{
				L_TRACE << "|";
			}
		}

		--total_reads_read;

		L_TRACE << this->results_to_string(total_reads_read);
	}

	string TagsFinderBase::results_to_string(long total_reads_read) const
	{
		stringstream ss;
		ss << " (" << total_reads_read << " reads)\n"
		   << this->get_additional_stat(total_reads_read) << "\n"
		   << this->trims_counter.print();

		return ss.str();
	}

	void TagsFinderBase::trim(const string &barcodes_tail, string &sequence, string &quality)
	{
		if (sequence.length() != quality.length())
			throw std::runtime_error("Read has different lengths of sequence and quality string: '" +
											 sequence + "', '" + quality + "'");

		len_t trim_pos = sequence.length();
		// attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
		// RC of UMI+second barcode (up to a length r1_rc_length - spacer_finder parameter)
		string rcb = Tools::reverse_complement(barcodes_tail);

		L_DEBUG << "-- barcode RC: " << rcb;
		len_t rc_pos = sequence.find(rcb);
		if (rc_pos != string::npos)
		{
			trim_pos = rc_pos;
			this->trims_counter.inc(TrimsCounter::RC);
			L_DEBUG << "-- found barcode RC at " << rc_pos;
		}
		else
		{
			// attempt 2: find polyA block
			rc_pos = sequence.find(this->poly_a);
			if (rc_pos != string::npos)
			{
				trim_pos = rc_pos;
				this->trims_counter.inc(TrimsCounter::POLY_A);
				L_DEBUG << "-- found polyA at " << rc_pos;
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
			this->trims_counter.inc(TrimsCounter::A_TRIM);
		}

		L_DEBUG << string(skip_count, '-') << "   trimming " << (sequence.length() - trim_pos);

		//attempt 4: apply
		if (sequence.length() != trim_pos)
		{
			sequence = sequence.substr(0, trim_pos);
			quality = quality.substr(0, trim_pos);
		}
		else
		{
			this->trims_counter.inc(TrimsCounter::NO_TRIM);
		}
	}
}
