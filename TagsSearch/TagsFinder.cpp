#include "TagsFinder.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "FilesProcessor.h"
#include "SpacerFinder.h"
#include "Tools/Logs.h"
#include "Tools/ReadsParameters.h"

using namespace std;

namespace TagsSearch
{
	TagsFinder::TagsFinder(const SpacerFinder &spacer_finder, const boost::property_tree::ptree &config)
			: spacer_finder(spacer_finder), trims_counter()
	{
		this->max_reads = config.get<long>("max_reads", -1);
		this->min_align_len = config.get<unsigned>("min_align_length");
		this->poly_a = config.get<string>("poly_a_tail");

		if (this->max_reads < 0)
		{
			this->max_reads = numeric_limits<long>::max();
		}
	}

	string TagsFinder::reverse_complement(const string &s)
	{
		char rcs[s.length()];

		for (int i = 0; i < s.length(); i++)
		{
			switch (s.at(s.length() - i - 1))
			{
				case 'A':
					rcs[i] = 'T';
					break;
				case 'T':
					rcs[i] = 'A';
					break;
				case 'C':
					rcs[i] = 'G';
					break;
				case 'G':
					rcs[i] = 'C';
					break;
				default:
					rcs[i] = 'N';
					break;
			}
		}
		string rcss = string(rcs, s.length());
		return (rcss);
	}

	bool TagsFinder::read_blocks(FilesProcessor &files_processor, long total_reads_read, string &out_1_line_2,
								 string &out_2_line_1, string &out_2_line_2, string &out_2_line_3, string &out_2_line_4)
	{
		string tmp;
		if (!files_processor.get_r1_line(tmp))
		{
			return 0;
		}

		if (tmp.at(0) != '@')
		{
			L_ERR << "read " << total_reads_read << ": R1 fastq malformed!";
			return 0;
		}

		if (!files_processor.get_r1_line(out_1_line_2))
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

		if (!files_processor.get_r2_line(out_2_line_1))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}

		if (out_2_line_1.at(0) != '@')
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq malformed!";
			return 0;
		}

		if (!files_processor.get_r2_line(out_2_line_2))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}
		if (!files_processor.get_r2_line(out_2_line_3))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}
		if (!files_processor.get_r2_line(out_2_line_4))
		{
			L_ERR << "read " << total_reads_read << ": R2 fastq ended prematurely!";
			return 0;
		}

		return 1;
	}

	void TagsFinder::run(const std::string &r1_filename, const std::string &r2_filename, const std::string &base_name,
						 bool save_reads_names)
	{
		L_TRACE << "reading reads ";

		FilesProcessor files_processor(r1_filename, r2_filename, base_name, this->max_reads);

		long total_reads_read = 1;
		Tools::reads_params_map_t reads_params;

		string r1_line2, r2_line1, r2_line2, r2_line3, r2_line4;
		while (TagsFinder::read_blocks(files_processor, total_reads_read, r1_line2, r2_line1, r2_line2, r2_line3, r2_line4))
		{
			if (total_reads_read % 1000000 == 0)
			{
				L_TRACE << "Total " << total_reads_read << " read";
			}

			Tools::ReadsParameters params = this->process_lines(total_reads_read++, r1_line2, r2_line2, r2_line3, r2_line4);
			if (params.is_empty)
				continue;

			string text;
			if (save_reads_names)
			{
				text = r2_line1;
//				reads_params[r2_line1] = params;
				reads_params.emplace(r2_line1, params);
			}
			else
			{
				text = params.to_string();
			}

			text += "\n" + r2_line2 + "\n" + r2_line3 + "\n" + r2_line4 + "\n";

			bool new_file = files_processor.write(text);

			if (new_file)
			{
				L_TRACE << "|";
			}
		}

		files_processor.close();
		if (save_reads_names)
		{
			files_processor.serialize_reads_params(reads_params);
		}

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

	TagsFinder::len_t TagsFinder::get_trim_position(len_t spacer_pos, const string &r1_line, const string &r2_line)
	{
		len_t r2_trim = r2_line.length();
		// attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
		// RC of UMI+second barcode (up to a length r1_rc_length - spacer_finder parameter)
		string rcb = this->spacer_finder.parse_r1_rc(r1_line, spacer_pos);
		rcb = TagsFinder::reverse_complement(rcb);

		L_DEBUG << "-- barcode RC: " << rcb;
		len_t rc_pos = r2_line.find(rcb);
		if (rc_pos != string::npos)
		{
			r2_trim = rc_pos;
			this->trims_counter.inc(TrimsCounter::RC);
			L_DEBUG << "-- found barcode RC at " << rc_pos;
		}
		else
		{
			// attempt 2: find polyA block
			rc_pos = r2_line.find(this->poly_a);
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
		while (r2_trim > 0 && (r2_line.at(r2_trim - 1) == 'A' || r2_line.at(r2_trim - 1) == 'N'))
		{
			r2_trim--;
			skip_count++;
			a_trim = true;
		}
		if (a_trim)
		{
			this->trims_counter.inc(TrimsCounter::A_TRIM);
		}

		L_DEBUG << string(skip_count, '-') << "   trimming " << (r2_line.length() - r2_trim);

		return r2_trim;
	}

	Tools::ReadsParameters TagsFinder::process_lines(long total_reads_read, const string &r1_line2, string &r2_line2,
									 const string &r2_line3, string &r2_line4)
	{
		L_DEBUG << r1_line2 << ":";

		SpacerFinder::len_t spacer_pos = this->spacer_finder.find_spacer(r1_line2);
		if (spacer_pos == SpacerFinder::ERR_CODE)
			return Tools::ReadsParameters();

		string cell_barcode = this->spacer_finder.parse_cell_barcode(r1_line2, spacer_pos);
		L_DEBUG << "-- cell barcode: " << cell_barcode << " (" << cell_barcode.length() << "nt)";

		string umi_barcode = this->spacer_finder.parse_umi_barcode(r1_line2, spacer_pos);

		L_DEBUG << "-- umi barcode: " << umi_barcode;
		L_DEBUG << "R2: " << r2_line2;

		// clean up R2
		len_t r2_trim = this->get_trim_position(spacer_pos, r1_line2, r2_line2);

		if (r2_line2.length() != r2_trim)
		{
			r2_line2 = r2_line2.substr(0, r2_trim);
			r2_line4 = r2_line4.substr(0, r2_trim);
		}
		else
		{
			this->trims_counter.inc(TrimsCounter::NO_TRIM);
		}

		L_DEBUG << " trimmed:" << r2_line2;

		// output
		if (r2_trim > this->min_align_len)
		{
			ostringstream text;
			text << '@' << total_reads_read << '!' << cell_barcode << '#' << umi_barcode << "\n";
			text << r2_line2 << "\n" << r2_line3 << "\n" << r2_line4 << endl;
			return Tools::ReadsParameters(total_reads_read, cell_barcode, umi_barcode);
		}

		return Tools::ReadsParameters();
	}

}