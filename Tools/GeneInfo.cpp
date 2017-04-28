#include "GeneInfo.h"

namespace Tools
{
	GeneInfo::GeneInfo(const std::string &chr_name, const std::string &id, const std::string &name, coord_t start_pos, coord_t end_pos)
			: Interval(start_pos, end_pos)
			, _chr_name(chr_name)
			, _id(id)
			, _name(name == id ? "" : name)
	{}

	bool GeneInfo::is_valid() const
	{
		return this->id() != "";
	}

	bool GeneInfo::operator<(const GeneInfo &other) const
	{
		return this->id() < other.id();
	}

	const std::string& GeneInfo::chr_name() const
	{
		return this->_chr_name;
	}

	const std::string& GeneInfo::id() const
	{
		return this->_id;
	}

	const std::string& GeneInfo::name() const
	{
		return this->_name == "" ? this->_id : this->_name;
	}

	GeneInfo::GeneInfo()
		: Interval(0, 0)
	{}
}