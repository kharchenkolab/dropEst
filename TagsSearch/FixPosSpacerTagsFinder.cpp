#include "FixPosSpacerTagsFinder.h"

#include <Tools/UtilFunctions.h>

#include <boost/algorithm/string.hpp>
#include <numeric>

namespace TagsSearch
{
	FixPosSpacerTagsFinder::FixPosSpacerTagsFinder(const std::string &barcode_fastq_name, const std::string &gene_fastq_name,
												   const boost::property_tree::ptree &barcodes_config,
												   const boost::property_tree::ptree &trimming_config, bool save_stats)
		: TagsFinderBase(trimming_config, save_stats)
		, _mask_parts(FixPosSpacerTagsFinder::parse_mask(barcodes_config.get<std::string>("barcode_mask", ""),
														 barcodes_config.get<std::string>("spacer_edit_dists", "")))
		, _trim_tail_length(std::min(barcodes_config.get<size_t>("r1_rc_length"),
									 std::accumulate(this->_mask_parts.begin(), this->_mask_parts.end(), (size_t)0,
													 [](size_t sum, const MaskPart & m) { return sum + m.length;})))
		, _outcomes(std::count_if(this->_mask_parts.begin(), this->_mask_parts.end(), [](const MaskPart & m){ return (m.type == MaskPart::SPACER);}))
		, _barcode_reader(barcode_fastq_name)
		, _gene_reader(gene_fastq_name)
	{}

	FixPosSpacerTagsFinder::MaskPart::MaskPart(const std::string &spacer, size_t length, Type type, size_t min_edit_distance)
			: spacer(spacer)
			, length(length)
			, type(type)
			, min_edit_distance(min_edit_distance)
	{}

	std::vector<FixPosSpacerTagsFinder::MaskPart> FixPosSpacerTagsFinder::parse_mask(const std::string& barcode_mask,
																					 const std::string& edit_dist_str)
	{
		std::vector<MaskPart> mask_parts;

		std::string mask = barcode_mask;
		std::string edit_dist_trim = edit_dist_str;

		boost::trim_if(mask, boost::is_any_of(" \t"));
		boost::trim_if(edit_dist_trim, boost::is_any_of(" \t"));
		if (mask.empty())
			throw std::runtime_error("Empty mask!");

		if (edit_dist_trim.empty())
			throw std::runtime_error("Empty edit distances!");

		std::vector<std::string> edit_distances;
		boost::split(edit_distances, edit_dist_trim, boost::is_any_of(", "), boost::token_compress_on);

		std::string::size_type cur_pos = 0;
		size_t spacer_ind = 0;
		while (cur_pos != mask.length())
		{
			std::string::size_type next_pos = mask.find_first_of("[(", cur_pos);
			if (next_pos == std::string::npos)
				throw std::runtime_error("Wrong mask format: " + mask);

			if (next_pos != cur_pos)
			{
				if (spacer_ind == edit_distances.size())
					throw std::runtime_error("Number of edit distances must be equal to the number of spacers");

				size_t ed = std::strtoul(edit_distances[spacer_ind++].c_str(), nullptr, 10);
				mask_parts.push_back(MaskPart(mask.substr(cur_pos, next_pos - cur_pos), next_pos - cur_pos, MaskPart::Type::SPACER, ed));
				cur_pos = next_pos;
			}
			++cur_pos;

			mask_parts.push_back(MaskPart());
			cur_pos = FixPosSpacerTagsFinder::parse_barcode_mask(mask, cur_pos, mask_parts.back()) + 1;
		}

		return mask_parts;
	}

	size_t FixPosSpacerTagsFinder::parse_barcode_mask(const std::string &mask, size_t cur_pos, MaskPart &mask_part)
	{
		FixPosSpacerTagsFinder::MaskPart::Type part_type;
		char end_char;

		if (mask[cur_pos - 1] == '[')
		{
			end_char = ']';
			part_type = MaskPart::CB;
		}
		else
		{
			end_char = ')';
			part_type = MaskPart::UMI;
		}
		size_t next_pos = mask.find(end_char, cur_pos);

		if (next_pos == std::string::npos)
			throw std::runtime_error("Wrong mask format: " + mask);

		mask_part.length = std::strtoul(mask.substr(cur_pos, next_pos - cur_pos).c_str(), nullptr, 10);
		mask_part.type = part_type;
		return next_pos;
	}

	bool FixPosSpacerTagsFinder::parse_fastq_record(FastQReader::FastQRecord &gene_record, Tools::ReadParameters &read_params)
	{
		auto barcodes_record = this->_barcode_reader.get_next_record();
		if (barcodes_record.id.empty())
			return false;

		gene_record = this->_gene_reader.get_next_record();
		if (gene_record.id.empty())
			throw std::runtime_error("File '" + this->_gene_reader.filename() + "', read '" + gene_record.id + "': fastq ended prematurely!");

		size_t seq_end = this->parse(barcodes_record.sequence, barcodes_record.quality, read_params);
		if (seq_end == std::string::npos)
		{
			read_params = Tools::ReadParameters();
			return true;
		}

		if (this->_trim_tail_length != 0)
		{
			auto tail = barcodes_record.sequence.substr(seq_end - this->_trim_tail_length, this->_trim_tail_length);
			this->trim(tail, gene_record.sequence, gene_record.quality);
		}

		return true;
	}

	size_t FixPosSpacerTagsFinder::parse(const std::string &r1_seq, const std::string &r1_quality, Tools::ReadParameters &read_params)
	{
		size_t cur_pos = 0, spacer_ind = 0;
		std::string cb, umi, cb_quality, umi_quality;
		for (auto const &mask_part : this->_mask_parts)
		{
			if (cur_pos + mask_part.length > r1_seq.length())
			{
				this->_outcomes.inc(MultiSpacerOutcomesCounter::SHORT_SEQ);
				return std::string::npos;
			}

			switch (mask_part.type)
			{
				case MaskPart::CB:
					cb += r1_seq.substr(cur_pos, mask_part.length);
					cb_quality += r1_quality.substr(cur_pos, mask_part.length);
					break;
				case MaskPart::SPACER:
					if (Tools::edit_distance(mask_part.spacer.c_str(), r1_seq.substr(cur_pos, mask_part.length).c_str(),
					                         mask_part.min_edit_distance) > mask_part.min_edit_distance)
					{
						this->_outcomes.inc_no_spacer(spacer_ind);
						return std::string::npos;
					}
					++spacer_ind;
					break;
				case MaskPart::UMI:
					umi += r1_seq.substr(cur_pos, mask_part.length);
					umi_quality += r1_quality.substr(cur_pos, mask_part.length);
					break;
				default:
					throw std::runtime_error("Unexpected MaskPart type: " + std::to_string(mask_part.type));
			}
			cur_pos += mask_part.length;
		}

		this->_outcomes.inc(MultiSpacerOutcomesCounter::OK);
		read_params = Tools::ReadParameters(cb, umi, cb_quality, umi_quality);
		return cur_pos;
	}

	std::string FixPosSpacerTagsFinder::get_additional_stat(long total_reads_read) const
	{
		return this->_outcomes.print(total_reads_read);
	}
}
