#include "TagsFinder.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "FilesProcessor.h"
#include "SpacerFinder.h"
#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"
#include "Tools/UtilFunctions.h"

using namespace std;

namespace TagsSearch
{
	TagsFinder::TagsFinder(const SpacerFinder &spacer_finder, const boost::property_tree::ptree &config)
			: spacer_finder(spacer_finder), trims_counter()
	{
		srand(time(nullptr));
		this->max_reads = config.get<long>("max_reads", -1);
		this->min_align_len = config.get<unsigned>("min_align_length");
		this->poly_a = config.get<string>("poly_a_tail");

		if (this->max_reads < 0)
		{
			this->max_reads = numeric_limits<long>::max();
		}
	}

	bool TagsFinder::read_blocks(FilesProcessor &files_processor, long total_reads_read, string &r1_seq,
								 string &r2_id, string &r2_seq, string &r2_description, string &r2_quality)
	{
		string tmp;
		if (!files_processor.get_r1_line(tmp))
			return 0;

		if (tmp.at(0) != '@')
		{
			L_ERR << "read " << total_reads_read << ": R1 fastq malformed!";
			return 0;
		}

		if (!files_processor.get_r1_line(r1_seq))
		{
			L_ERR << "read " << total_reads_read << ": R1 fastq ended prematurely!";
			return 0;
		}
		if (!files_processor.get_r1_line(tmp))
		{
			L_ERR << "read " << total_reads_read << ": R1 fastq ended prematurely!";
			return 0;
		}
		if (!files_processor.get_r1_line(tmp))
		{
			L_ERR << "read " << total_reads_read << ": R1 fastq ended prematurely!";
			return 0;
		}

		if (!files_processor.get_r2_line(r2_id))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}

		if (r2_id.at(0) != '@')
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq malformed!";
			return 0;
		}

		if (!files_processor.get_r2_line(r2_seq))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}
		if (!files_processor.get_r2_line(r2_description))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}
		if (!files_processor.get_r2_line(r2_quality))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}
		if (r2_seq.length() != r2_quality.length())
		{
			L_ERR << "read " << total_reads_read << " have qual length (" << r2_quality.length() <<
						") differs from seq length (" << r2_seq.length() << ")";
			return 0;
		}

		return 1;
	}

	void TagsFinder::run(const std::string &r1_filename, const std::string &r2_filename, const std::string &base_name,
						 bool save_reads_names)
	{
		L_TRACE << "reading reads ";

		FilesProcessor files_processor(r1_filename, r2_filename, base_name, this->max_reads);

		long total_reads_read = 1, parsed_reads = 0;
		std::string file_uid = std::to_string(rand()) + char(rand() % 25 + 'A');

		string r1_seq, r2_id, r2_seq, r2_description, r2_quality_str;
		while (TagsFinder::read_blocks(files_processor, total_reads_read, r1_seq, r2_id, r2_seq, r2_description,
									   r2_quality_str))
		{
			if (total_reads_read % 1000000 == 0)
			{
				L_TRACE << "Total " << total_reads_read << " read (" << parsed_reads << " parsed)";
			}

			Tools::ReadParameters params = this->parse_and_trim(r1_seq, r2_id, r2_seq, r2_quality_str);

			++total_reads_read;
			if (params.is_empty())
				continue;

			++parsed_reads;

			string text;
			std::string new_id = "@" + file_uid + std::to_string(total_reads_read);
			if (save_reads_names)
			{
				text = new_id;
				files_processor.write_read_params(new_id, params);
			}
			else
			{
				text = params.to_monolithic_string(new_id);
			}

			text += "\n" + r2_seq + "\n" + r2_description + "\n" + r2_quality_str + "\n";

			bool new_file = files_processor.write(text);

			if (new_file)
			{
				L_TRACE << "|";
			}
		}

		files_processor.close();

		--total_reads_read;

		L_TRACE << this->results_to_string(total_reads_read);
	}

	string TagsFinder::results_to_string(long total_reads_read) const
	{
		stringstream ss;
		ss << " (" << total_reads_read << " reads)";
		ss << this->spacer_finder.get_outcomes_counter().print(total_reads_read);
		ss << this->trims_counter.print();

		return ss.str();
	}

	TagsFinder::len_t TagsFinder::get_trim_position(len_t spacer_end, const string &r1_seq, const string &r2_seq)
	{
		len_t r2_trim = r2_seq.length();
		// attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
		// RC of UMI+second barcode (up to a length r1_rc_length - spacer_finder parameter)
		string rcb = this->spacer_finder.parse_r1_rc(r1_seq, spacer_end);
		rcb = Tools::reverse_complement(rcb);

		L_DEBUG << "-- barcode RC: " << rcb;
		len_t rc_pos = r2_seq.find(rcb);
		if (rc_pos != string::npos)
		{
			r2_trim = rc_pos;
			this->trims_counter.inc(TrimsCounter::RC);
			L_DEBUG << "-- found barcode RC at " << rc_pos;
		}
		else
		{
			// attempt 2: find polyA block
			rc_pos = r2_seq.find(this->poly_a);
			if (rc_pos != -1)
			{
				r2_trim = rc_pos;
				this->trims_counter.inc(TrimsCounter::POLY_A);
				L_DEBUG << "-- found polyA at " << rc_pos;
			}
		}

		// attempt 3: trim trailing As
		bool a_trim = false;
		len_t skip_count = 0;
		while (r2_trim > 0 && (r2_seq.at(r2_trim - 1) == 'A' || r2_seq.at(r2_trim - 1) == 'N'))
		{
			r2_trim--;
			skip_count++;
			a_trim = true;
		}
		if (a_trim)
		{
			this->trims_counter.inc(TrimsCounter::A_TRIM);
		}

		L_DEBUG << string(skip_count, '-') << "   trimming " << (r2_seq.length() - r2_trim);

		return r2_trim;
	}

	Tools::ReadParameters TagsFinder::parse_and_trim(const string &r1_seq, const string &r2_id, string &r2_seq,
													 string &r2_quality_str)
	{
		if (r2_seq.length() != r2_quality_str.length())
			throw std::runtime_error("Read " + r2_id + " have different length of sequence and quality string: '" +
											 r2_seq + "', '" + r2_quality_str + "'");

		L_DEBUG << r1_seq << ":";

		auto spacer_pos = this->spacer_finder.find_spacer(r1_seq);
		if (spacer_pos.first == SpacerFinder::ERR_CODE)
			return Tools::ReadParameters();

		string cell_barcode = this->spacer_finder.parse_cell_barcode(r1_seq, spacer_pos.first, spacer_pos.second);
		L_DEBUG << "-- cell barcode: " << cell_barcode << " (" << cell_barcode.length() << "nt)";

		string umi_barcode = this->spacer_finder.parse_umi_barcode(r1_seq, spacer_pos.second);
		if (umi_barcode.length() == 0)
		{
			umi_barcode = "N";
		}

		L_DEBUG << "-- umi barcode: " << umi_barcode;
		L_DEBUG << "R2: " << r2_seq;

		len_t r2_trim = this->get_trim_position(spacer_pos.second, r1_seq, r2_seq);

		if (r2_seq.length() != r2_trim)
		{
			r2_seq = r2_seq.substr(0, r2_trim);
			r2_quality_str = r2_quality_str.substr(0, r2_trim);
		}
		else
		{
			this->trims_counter.inc(TrimsCounter::NO_TRIM);
		}

		L_DEBUG << " trimmed:" << r2_seq;

		if (r2_trim < this->min_align_len)
			return Tools::ReadParameters();

		return Tools::ReadParameters(r2_id, cell_barcode, umi_barcode);
	}

}
