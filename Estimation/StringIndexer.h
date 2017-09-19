#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace Estimation
{
	class StringIndexer
	{
	public:
		typedef unsigned index_t;
		typedef std::vector<std::string> values_t;
	private:
		values_t _values;
		std::unordered_map<std::string, index_t> _indexes;

	public:
		const values_t& values() const;
		const std::string& get_value(index_t index) const;
		index_t get_index(const std::string &value) const;
		index_t add(const std::string& value);
	};
}