#include "MergeUMIsStrategyDirectional.h"

#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>

namespace Estimation
{
namespace Merge
{
namespace UMIs
{

	MergeUMIsStrategyDirectional::MergeUMIsStrategyDirectional(double mult, unsigned max_edit_distance)
		: _mult(mult)
		, _max_edit_distance(max_edit_distance)
	{}

	void MergeUMIsStrategyDirectional::merge(CellsDataContainer &container) const
	{
		Tools::trace_time("Start UMI merge");
		size_t total_cell_merged = 0, total_umi_merged = 0, total_cells_processed = 0;
		for (size_t cell_id = 0; cell_id < container.total_cells_number(); ++cell_id)
		{
			auto const &cell = container.cell(cell_id);
			if (!cell.is_real())
				continue;

			total_cells_processed++;
			for (auto const &gene : cell.genes())
			{
				umi_vec_t umis;
				for (auto const &umi : gene.second.umis())
				{
					auto umi_seq = container.umi_indexer().get_value(umi.first);
					umis.emplace_back(umi_seq, umi.second.read_count());
				}


				auto merge_targets = this->find_targets(umis);
				if (merge_targets.empty())
					continue;

				total_cell_merged++;
				total_umi_merged += merge_targets.size();

				container.merge_umis(cell_id, gene.first, merge_targets);
			}
		}

		L_TRACE << total_cells_processed << " cells processed. Merged " << total_umi_merged << " UMIs from "
		        << total_cell_merged << " cells.";
		Tools::trace_time("UMI merge finished");
	}

	MergeUMIsStrategyDirectional::merge_targets_t MergeUMIsStrategyDirectional::find_targets(umi_vec_t &umis) const
	{
		std::sort(umis.begin(), umis.end(), [](const UmiWrap &u1, const UmiWrap &u2){return u1.n_reads < u2.n_reads;});
		CellsDataContainer::s_s_hash_t merge_targets;

		for (size_t src_id = 0; src_id < umis.size(); ++src_id)
		{
			std::string target = this->find_target(src_id, umis);
			if (!target.empty())
			{
				merge_targets[umis[src_id].sequence] = target;
			}
		}

		for (long i = umis.size() - 1; i >= 0; --i)
		{
			auto dst_iter = merge_targets.find(umis[i].sequence);
			if (dst_iter == merge_targets.end())
				continue;

			dst_iter = merge_targets.find(dst_iter->second);
			if (dst_iter == merge_targets.end())
				continue;

			merge_targets[umis[i].sequence] = dst_iter->second;
		}

		return merge_targets;
	}

	std::string MergeUMIsStrategyDirectional::find_target(size_t src_id, MergeUMIsStrategyDirectional::umi_vec_t &umis) const
	{
		auto const &src_umi = umis[src_id];
		const bool has_ns = (src_umi.sequence.find('N') != std::string::npos);

		std::string target;
		unsigned min_ed = std::numeric_limits<unsigned>::max();
		for (long dst_id = umis.size() - 1; dst_id > src_id; --dst_id)
		{
			auto const &dst_umi = umis[dst_id];
			if (src_umi.n_reads * this->_mult > dst_umi.n_reads)
				break;

			auto ed = Tools::edit_distance(src_umi.sequence.c_str(), dst_umi.sequence.c_str(), true, this->_max_edit_distance);
			if (ed > this->_max_edit_distance)
				continue;

			if (ed < min_ed)
			{
				target = dst_umi.sequence;
				if (!has_ns && ed <= 1 || ed == 0)
					break;

				min_ed = ed;
			}
		}

		if (has_ns && target.empty())
			return MergeUMIsStrategyAbstract::fix_n_umi_with_random(src_umi.sequence);

		return target;
	}

	MergeUMIsStrategyDirectional::UmiWrap::UmiWrap(const std::string &sequence, size_t n_reads)
		: sequence(sequence)
		, n_reads(n_reads)
	{}
}
}
}