#include "GeneInfo.h"

namespace Tools
{
	GeneInfo::GeneInfo(const std::string &chr_name, std::string id, pos_t start_pos, pos_t end_pos)
			: _chr_name(chr_name)
			, _id(id)
			, _start_pos(start_pos)
			, _end_pos(end_pos)
	{
		this->_chr_num = GeneInfo::parse_chr_name(chr_name);
	}

	bool GeneInfo::is_valid() const
	{
		return this->id() != "";
	}

	bool GeneInfo::operator<(const GeneInfo &other) const
	{
		return this->id() < other.id();
	}

	bool GeneInfo::is_intercept(const GeneInfo &other) const
	{
		return this->start_pos() <= other.end_pos() && this->end_pos() >= other.start_pos();
	}

	void GeneInfo::merge(const GeneInfo &other)
	{
		this->_start_pos = std::min(this->_start_pos, other._start_pos);
		this->_end_pos= std::max(this->_end_pos, other._end_pos);
	}

	std::string GeneInfo::chr_name() const
	{
		return this->_chr_name;
	}

	std::string GeneInfo::id() const
	{
		return this->_id;
	}

	GeneInfo::pos_t GeneInfo::start_pos() const
	{
		return this->_start_pos;
	}

	GeneInfo::pos_t GeneInfo::end_pos() const
	{
		return this->_end_pos;
	}

	GeneInfo::pos_t GeneInfo::size() const
	{
		return this->end_pos() - this->start_pos();
	}

	GeneInfo::num_t GeneInfo::chr_num() const
	{
		return this->_chr_num;
	}

	GeneInfo::num_t GeneInfo::parse_chr_name(const std::string &chr_name) //TODO Replace index by full name
	{
		return chr_name.length() > 2 ? (num_t)atoi(chr_name.substr(3).c_str()) : (num_t)atoi(chr_name.c_str());
	}
}