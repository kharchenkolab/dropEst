#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace Estimation
{
	class StringIndexer
	{
	public:
		using index_t=unsigned;
	private:
		std::vector<std::string> _values;
		std::unordered_map<std::string, index_t> _indexes;

	public:
		const std::string& get_value(index_t index) const;
		index_t get_index(const std::string &value) const;
		index_t add(const std::string& value);
	};
}