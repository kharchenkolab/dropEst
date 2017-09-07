#include "StringIndexer.h"

namespace Estimation
{
	const std::string &StringIndexer::get_value(index_t index) const
	{
		return this->_values.at(index);
	}

	StringIndexer::index_t StringIndexer::add(const std::string &value)
	{
		auto iter = this->_indexes.emplace(value, this->_indexes.size());
		if (iter.second)
		{
			this->_values.push_back(value);
		}
		return iter.first->second;
	}

	StringIndexer::index_t StringIndexer::get_index(const std::string &value) const
	{
		return this->_indexes.at(value);
	}

}
