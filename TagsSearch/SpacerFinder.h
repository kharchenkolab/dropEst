#pragma once

#include <string>

#include <boost/property_tree/ptree.hpp>

#include "Counters/OutcomesCounter.h"

namespace Tools
{
	class ReadParameters;
}

namespace TagsSearch
{
// inDrop read structure: [C_(spacer_min_pos-spacer_max_pos)][Spacer][C(barcode_length)][M(umi_length)]
	class SpacerFinder
	{
	public:
		typedef std::string::size_type len_t;
		typedef std::pair<len_t, len_t> spacer_pos_t;
		static const len_t ERR_CODE = std::string::npos;

	private:
		std::string spacer_prefix;
		std::string spacer_suffix;

		std::string spacer;
		size_t max_spacer_ed;
		size_t spacer_min_pos;
		size_t spacer_max_pos;
		len_t barcode_length;
		len_t umi_length;
		len_t r1_rc_length;

		size_t spacer_max_suffix_start;
		size_t spacer_min_suffix_start;

		size_t min_seq_len;

		OutcomesCounter outcomes;

	public:
		SpacerFinder()
		{}

		SpacerFinder(const boost::property_tree::ptree &config, const std::string& reads_params_file = "");

		spacer_pos_t find_spacer(const std::string& seq);

		std::string parse_cell_barcode(const std::string& seq, len_t spacer_start, len_t spacer_end) const;
		std::string parse_umi_barcode(const std::string& seq, len_t spacer_end) const;
		std::string parse_r1_rc(const std::string &seq, len_t spacer_end) const;

		const OutcomesCounter& get_outcomes_counter() const;

	private:
		spacer_pos_t find_spacer_partial(const std::string& seq);
	};
}