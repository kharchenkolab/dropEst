#pragma once

#include <string>
#include "Interval.h"

namespace TestTools
{
	struct testGtf;
	struct testGeneMerge;
}

namespace Tools
{
	class GeneInfo : public Interval
	{
		friend struct TestTools::testGtf;
		friend struct TestTools::testGeneMerge;

	private:
		std::string _chr_name = "";
		std::string _id = "";
		std::string _name = "";

	public:
		GeneInfo();
		GeneInfo(const std::string &chr_name, const std::string &id, const std::string &name, coord_t start_pos, coord_t end_pos);

		const std::string& chr_name() const;
		const std::string& id() const;
		const std::string& name() const;

		bool is_valid() const;
		bool operator<(const GeneInfo &other) const;
	};
}