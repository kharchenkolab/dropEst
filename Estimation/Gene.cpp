#include "Gene.h"

namespace Estimation
{
	const UMI &Estimation::Gene::at(const std::string &umi) const
	{
		return this->_umis.at(umi);
	}

	bool Gene::add_umi(const std::string &umi, const UMI::Mark &umi_mark)
	{
		auto insert_it = this->_umis.emplace(umi, 0);
		auto &new_umi = insert_it.first->second;
		new_umi.read_count++;
		new_umi.mark.add(umi_mark);

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

		auto source_umi_it = this->_umis.find(source_umi);
		auto target_umi_it = this->_umis.emplace(target_umi, source_umi_it->second);
		if (!target_umi_it.second)
		{
			target_umi_it.first->second.merge(this->_umis.at(source_umi));
		}
		this->_umis.erase(source_umi_it);
	}

	size_t Gene::number_of_requested_umis(const UMI::Mark::query_t &query, bool return_reads) const
	{
		size_t requested_umis_num = 0;
		for (auto const &umi : this->_umis)
		{
			if (!umi.second.mark.match(query))
				continue;

			if (return_reads)
			{
				requested_umis_num += umi.second.read_count;
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
			reads_num += umi.second.read_count;
		}

		return reads_num;
	}

	Gene::s_ul_hash_t Gene::requested_reads_per_umi(const UMI::Mark::query_t &query) const
	{
		s_ul_hash_t reads_per_umi;
		for (auto const &umi : this->_umis)
		{
			if (!umi.second.mark.match(query))
				continue;

			reads_per_umi.emplace(umi.first, umi.second.read_count);
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
		return this->_umis.find(umi) != this->_umis.end();
	}
}