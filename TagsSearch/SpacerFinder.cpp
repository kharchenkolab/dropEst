#include "SpacerFinder.h"

#include "Tools/UtilFunctions.h"
#include "Tools/Logs.h"

using std::string;

namespace TagsSearch
{
	SpacerFinder::SpacerFinder(const boost::property_tree::ptree &config, const std::string &reads_params_file)
			: reads_params_file(reads_params_file), need_save_reads_params(reads_params_file.length() != 0)
	{
		this->spacer = config.get<string>("spacer");
		this->spacer_2 = config.get<string>("spacer2");
		this->max_spacer_ed = config.get<unsigned>("max_spacer_edit_distance");
		this->spacer_prefix_length = config.get<unsigned>("spacer_prefix_length");

		this->spacer_min_pos = config.get<unsigned>("spacer_min_pos");
		this->spacer_max_pos = config.get<unsigned>("spacer_max_pos");
		this->barcode_length = config.get<unsigned>("barcode_length");
		this->umi_length = config.get<unsigned>("umi_length");
		this->r1_rc_length = config.get<unsigned>("r1_rc_length");

		this->min_seq_len = this->spacer_min_pos + this->barcode_length + this->umi_length + this->spacer.length();

		if (this->spacer.length() != this->spacer_2.length())
			throw std::runtime_error("Spacers must be the same length: '" + spacer + "' '" + spacer_2 + "'");

		if (this->spacer.length() <= SpacerFinder::spacer_prefix_length)
			throw std::runtime_error(
					"Spacers length should be bigger then spacer_prefix_length ('" + spacer + "): '" + spacer + "'");
	}

	string::size_type SpacerFinder::find_spacer(const string &seq)
	{
		if (seq.length() < this->min_seq_len)
		{
			this->outcomes.inc(OutcomesCounter::SHORT_SEQ);
			L_DEBUG << "-- read is too short";

			return string::npos;
		}

		len_t spacer_pos = this->find_spacer_full(seq);
		if (spacer_pos == string::npos)
		{
			spacer_pos = this->find_spacer_partial(seq);
		}

		if (spacer_pos == string::npos)
			return spacer_pos;

		if (spacer_pos < this->spacer_min_pos || spacer_pos > this->spacer_max_pos)
		{
			this->outcomes.inc(OutcomesCounter::SPACER_MISPLACED);
			L_DEBUG << "-- invalid spacer position";

			return string::npos;
		}

		if (seq.length() < spacer_pos + spacer.length() + this->barcode_length + this->umi_length + 1)
		{
			this->outcomes.inc(OutcomesCounter::SHORT_SEQ);
			L_DEBUG << "-- read is too short for this spacer";

			return string::npos;
		}

		this->outcomes.inc(OutcomesCounter::OK);
		return spacer_pos;
	}

	SpacerFinder::len_t SpacerFinder::find_spacer_partial(const string &seq)
	{
		len_t suffix_start = this->spacer.length() - this->spacer_prefix_length;
		len_t suffix_seq_pos = seq.rfind(this->spacer.substr(suffix_start, this->spacer_prefix_length));
		if (suffix_seq_pos != string::npos && suffix_seq_pos >= suffix_start)
		{
			len_t spacer_pos = suffix_seq_pos - suffix_start;

			if (spacer_pos >= this->spacer_min_pos && spacer_pos < this->spacer_max_pos)
			{
				int ed = Tools::edit_distance(this->spacer.c_str(),
											  seq.substr(spacer_pos, this->spacer.length()).c_str());

				L_DEBUG << "-- postfix match at " << suffix_seq_pos << " ed=" << ed;
				L_DEBUG << "-- given:" << this->spacer;
				L_DEBUG << "-- match:" << seq.substr(spacer_pos, this->spacer.length());

				if (ed <= SpacerFinder::max_spacer_ed)
				{
					this->outcomes.inc(OutcomesCounter::SPACER_2);
					return spacer_pos;
				}
			}
		}

		len_t prefix_pos = seq.rfind(this->spacer.substr(0, this->spacer_prefix_length), this->spacer_max_pos);

		if (prefix_pos != string::npos && prefix_pos < this->spacer_max_pos)
		{
			int ed = Tools::edit_distance(this->spacer.c_str(), seq.substr(prefix_pos, this->spacer.length()).c_str());
			L_DEBUG << "-- prefix match at " << prefix_pos << " ed=" << ed;

			if (ed <= SpacerFinder::max_spacer_ed)
			{
				this->outcomes.inc(OutcomesCounter::SPACER_2);
				return prefix_pos;
			}
		}

		L_DEBUG << "-- spacer not found";

		this->outcomes.inc(OutcomesCounter::NO_SPACER);
		return string::npos;
	}

	SpacerFinder::len_t SpacerFinder::find_spacer_full(const string &seq)
	{
		len_t spacer_pos = seq.find(this->spacer);
		if (spacer_pos != string::npos)
			return spacer_pos;

		spacer_pos = seq.find(this->spacer_2);
		if (spacer_pos != string::npos)
		{
			this->outcomes.inc(OutcomesCounter::SPACER_2);
			L_DEBUG << "-- secondary spacer" << std::endl;

			return spacer_pos;
		}
		return string::npos;
	}

	string SpacerFinder::parse_cell_barcode(const string &seq, len_t spacer_pos) const
	{
		return seq.substr(0, spacer_pos) + seq.substr(spacer_pos + this->spacer.length(), this->barcode_length);
	}

	string SpacerFinder::parse_umi_barcode(const string &seq, len_t spacer_pos) const
	{
		return seq.substr(spacer_pos + spacer.length() + this->barcode_length, this->umi_length);
	}

	string SpacerFinder::parse_r1_rc(const string &seq, len_t spacer_pos) const
	{
		return seq.substr(
				spacer_pos + this->spacer.length() + this->barcode_length + this->umi_length - this->r1_rc_length,
				this->r1_rc_length);
	}

	const OutcomesCounter &SpacerFinder::get_outcomes_counter() const
	{
		return this->outcomes;
	}

}