#pragma once

#include <string>

#include <boost/property_tree/ptree.hpp>

#include "Counters/OutcomesCounter.h"

// inDrop read structure: [C_(spacer_min_pos-spacer_max_pos)][Spacer][C(barcode_length)][M(umi_length)]
class SpacerFinder
{
public:
	typedef std::string::size_type len_t;
	static const len_t ERR_CODE = std::string::npos;

private:
	unsigned max_spacer_ed; // maximum number of spacer mismatches
	unsigned spacer_prefix_length; // length of the spacer prefix/suffix to use as a seed
	std::string spacer;
	std::string spacer_2;

	int spacer_min_pos;
	int spacer_max_pos;
	len_t barcode_length;
	len_t umi_length;
	len_t r1_rc_length;

	size_t min_seq_len;

	OutcomesCounter outcomes;

public:
	SpacerFinder() = default;
	SpacerFinder(const boost::property_tree::ptree &config);

	len_t find_spacer(const std::string& seq);

	std::string parse_cell_barcode(const std::string& seq, len_t spacer_pos) const;
	std::string parse_umi_barcode(const std::string& seq, len_t spacer_pos) const;
	std::string parse_r1_rc(const std::string &seq, len_t spacer_pos) const;

	const OutcomesCounter& get_outcomes_counter() const;

private:
	len_t find_spacer_partial(const std::string& seq);
	len_t find_spacer_full(const std::string& seq);
};
