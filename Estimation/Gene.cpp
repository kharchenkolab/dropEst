#include "Gene.h"

#include "ReadInfo.h"

namespace Estimation
{
	Gene::Gene(StringIndexer *umi_indexer, bool save_merge_targets)
		: _umi_indexer(umi_indexer)
		, _save_merge_targets(save_merge_targets)
	{}

	const UMI &Estimation::Gene::at(const std::string &umi) const
	{
		return this->_umis.at(this->_umi_indexer->get_index(umi));
	}

	bool Gene::add_umi(const ReadInfo &read_info)
	{
		auto umi_index = this->_umi_indexer->add(read_info.params.umi());
		auto insert_it = this->_umis.emplace(umi_index, UMI(read_info.params.umi_quality().length(), 0));
		insert_it.first->second.add_read(read_info);

		return insert_it.second;
	}

	void Gene::merge(const Gene &source)
	{
		for (auto const &merged_umi: source._umis)
		{
			auto insert_it = this->_umis.insert(merged_umi);
			if (insert_it.second)
				continue;

			insert_it.first->second.merge(merged_umi.second);
		}
	}

	void Gene::merge(const std::string &source_umi, const std::string &target_umi)
	{
		if (source_umi == target_umi)
			return;

		auto source_umi_it = this->_umis.find(this->_umi_indexer->get_index(source_umi));
		if (source_umi_it == this->_umis.end())
			throw std::runtime_error("Source UMI doesn't belong to the gene: " + source_umi);

		auto target_umi_it = this->_umis.emplace(this->_umi_indexer->add(target_umi), source_umi_it->second);
		if (!target_umi_it.second)
		{
			target_umi_it.first->second.merge(source_umi_it->second);
		}
		this->_umis.erase(source_umi_it);

		if (this->_save_merge_targets)
		{
			this->_merge_targets[source_umi] = target_umi;
		}
	}

	size_t Gene::number_of_requested_umis(const UMI::Mark::query_t &query, bool return_reads) const
	{
		size_t requested_umis_num = 0;
		for (auto const &umi : this->_umis)
		{
			if (!umi.second.mark().match(query))
				continue;

			if (return_reads)
			{
				requested_umis_num += umi.second.read_count();
			}
			else
			{
				requested_umis_num++;
			}
		}

		return requested_umis_num;
	}

	size_t Gene::number_of_umis(bool return_reads) const
	{
		if (!return_reads)
			return this->_umis.size();

		size_t reads_num = 0;
		for (auto const &umi : this->_umis)
		{
			reads_num += umi.second.read_count();
		}

		return reads_num;
	}

	Gene::s_ul_hash_t Gene::requested_reads_per_umi(const UMI::Mark::query_t &query) const
	{
		s_ul_hash_t reads_per_umi;
		for (auto const &umi : this->_umis)
		{
			if (!umi.second.mark().match(query))
				continue;

			reads_per_umi.emplace(this->_umi_indexer->get_value(umi.first), umi.second.read_count());
		}

		return reads_per_umi;
	}

	const Gene::umis_t &Gene::umis() const
	{
		return this->_umis;
	}

	size_t Gene::size() const
	{
		return this->_umis.size();
	}

	bool Gene::has(const std::string &umi) const
	{
		return this->_umis.find(this->_umi_indexer->get_index(umi)) != this->_umis.end();
	}

	const Gene::s_s_hash_t &Gene::merge_targets() const
	{
		return this->_merge_targets;
	}
}