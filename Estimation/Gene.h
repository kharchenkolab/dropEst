#pragma once

#include <map>
#include <string>
#include <unordered_map>

#include "UMI.h"

namespace Estimation
{
	class Gene
	{
	private:
		typedef std::unordered_map<std::string, size_t> s_ul_hash_t;

	public:
		typedef std::map<std::string, UMI> umis_t;

	private:
		umis_t _umis;

	public:
		const UMI& at(const std::string& umi) const;
		const umis_t& umis() const;
		size_t size() const;
		bool has(const std::string& umi) const;

		size_t number_of_requested_umis(const UMI::Mark::query_t &query, bool return_reads) const;
		size_t number_of_umis(bool return_reads) const;
		s_ul_hash_t requested_reads_per_umi(const UMI::Mark::query_t &query) const;

		bool add_umi(const std::string &umi, const UMI::Mark &umi_mark);
		void merge(const Gene& source);
		void merge(const std::string& source_umi, const std::string& target_umi);
	};
}