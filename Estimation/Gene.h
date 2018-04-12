#pragma once

#include <map>
#include <string>
#include <unordered_map>

#include "UMI.h"
#include "StringIndexer.h"

namespace Estimation
{
	class ReadInfo;
	class Gene
	{
	private:
		using s_ul_hash_t = std::unordered_map<std::string, size_t>;

	public:
		using umis_t = std::map<StringIndexer::index_t, UMI>;
		using s_s_hash_t = std::unordered_map<std::string, std::string>;

	private:
		const bool _save_merge_targets;
		s_s_hash_t _merge_targets;

		umis_t _umis;
		StringIndexer *_umi_indexer;

	public:
		Gene(StringIndexer *umi_indexer, bool save_merge_targets);

		const UMI& at(const std::string& umi) const;
		const umis_t& umis() const;
		size_t size() const;
		bool has(const std::string& umi) const;
		const s_s_hash_t& merge_targets() const;

		size_t number_of_requested_umis(const UMI::Mark::query_t &query, bool return_reads) const;
		size_t number_of_umis(bool return_reads) const;
		s_ul_hash_t requested_reads_per_umi(const UMI::Mark::query_t &query) const;

		bool add_umi(const ReadInfo &read_info);
		void merge(const Gene& source);
		void merge(const std::string& source_umi, const std::string& target_umi);
	};
}
