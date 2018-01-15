#include <Tools/IndexedValue.h>
#include <Tools/UtilFunctions.h>
#include "BarcodesParser.h"

namespace Estimation
{
namespace Merge
{
namespace BarcodesParsing
{
	BarcodesParser::BarcodesParser(const std::string &barcodes_filename)
		: _barcodes_filename(barcodes_filename)
	{}

	BarcodesParser::BarcodesDistance::BarcodesDistance(const std::vector<size_t> &barcodes_inds,
													   unsigned int edit_distance)
			: barcode_part_inds(barcodes_inds)
			, edit_distance(edit_distance)
	{}

	BarcodesParser::edit_distance_parts_list_t BarcodesParser::get_distances_to_barcode(const std::string &barcode) const
	{
		barcodes_list_t barcode_parts = this->split_barcode(barcode);
		edit_distance_parts_list_t res(this->_barcodes.size());

		for (size_t part_ind = 0; part_ind < this->_barcodes.size(); ++part_ind)
		{
			auto &part_res = res[part_ind];
			auto &part_barcodes = this->_barcodes[part_ind];
			auto const &cb_part = barcode_parts[part_ind];
			for (size_t cb_ind = 0; cb_ind < this->_barcodes[part_ind].size(); cb_ind++)
			{
				part_res.emplace_back(cb_ind, Tools::edit_distance(cb_part.c_str(), part_barcodes[cb_ind].c_str()));
			}
			std::sort(part_res.begin(), part_res.end(), Tools::IndexedValue::value_less);
		}

		return res;
	}

	std::vector<BarcodesParser::BarcodesDistance> BarcodesParser::get_real_neighbour_cbs(const std::string &barcode) const
	{
		edit_distance_parts_list_t dist_parts = this->get_distances_to_barcode(barcode);
		std::vector<BarcodesDistance> res;
		if (dist_parts.empty())
			return res;

		this->push_remaining_dists(dist_parts.begin(), dist_parts.end(), 0, std::vector<size_t>(), res);
		return res;
	}

	void BarcodesParser::push_remaining_dists(edit_distance_parts_list_t::const_iterator begin,
											  edit_distance_parts_list_t::const_iterator end,
											  unsigned edit_distance, const std::vector<size_t> &barcodes_inds,
											  BarcodesParser::barcodes_distance_list_t &res) const
	{
		if (begin == end)
		{
			res.emplace_back(barcodes_inds, edit_distance);
			return;
		}

		std::vector<size_t> cur_inds(barcodes_inds);
		cur_inds.push_back(0);
		for (auto const &cur_dist : *begin)
		{
			unsigned cur_ed = edit_distance + cur_dist.value;
			if (cur_ed > BarcodesParser::MAX_REAL_MERGE_EDIT_DISTANCE)
				return;

			cur_inds.back() = cur_dist.index;
			this->push_remaining_dists(begin + 1, end, cur_ed, cur_inds, res);
		}
	}

	std::string BarcodesParser::get_barcode(const std::vector<size_t> &barcode_part_inds) const
	{
		if (barcode_part_inds.size() != this->_barcodes.size())
			throw std::runtime_error("Unexpected number of barcode parts: " + std::to_string(barcode_part_inds.size()));
		std::string res;
		for (size_t i = 0; i < barcode_part_inds.size(); ++i)
		{
			res += this->_barcodes[i].at(barcode_part_inds[i]);
		}
		return res;
	}

	void BarcodesParser::init()
	{
		this->_barcodes = this->get_barcodes_list(this->_barcodes_filename);

		if (this->_barcodes.empty())
			throw std::runtime_error("ERROR: empty barcodes list");

		for (auto const &barcode_parts : this->_barcodes)
		{
			if (barcode_parts.empty())
				throw std::runtime_error("ERROR: empty barcodes list");
		}
	}

	const std::string &BarcodesParser::barcode(size_t part_ind, size_t barcode_ind) const
	{
		return this->_barcodes.at(part_ind).at(barcode_ind);
	}

	void BarcodesParser::release()
	{
		this->_barcodes.clear();
	}

	const size_t BarcodesParser::barcode_parts_num() const
	{
		return this->_barcodes.size();
	}

	bool BarcodesParser::read_line(std::ifstream &barcodes_file, barcodes_list_t &barcodes, bool require_equal_length)
	{
		std::string line;
		Tools::ReverseComplement rc;
		if (!std::getline(barcodes_file, line))
			return false;

		std::istringstream b1_in(line);
		std::string::size_type barcode_length = 0;
		while (b1_in)
		{
			std::string barcode;
			b1_in >> barcode;
			if (barcode.empty())
				continue;

			if (barcode_length == 0)
			{
				barcode_length = barcode.length();
			}
			else if (require_equal_length && barcode_length != barcode.length())
				throw std::runtime_error("All barcodes in one line must have the same length");

			barcodes.push_back(rc.rc(barcode));
		}

		return true;
	}
}
}
}
