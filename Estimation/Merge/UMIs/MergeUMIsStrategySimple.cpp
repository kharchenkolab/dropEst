#include "MergeUMIsStrategySimple.h"

#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>

#include <limits>

namespace Estimation
{
namespace Merge
{
namespace UMIs
{
const std::string MergeUMIsStrategySimple::nucleotides = "ACGT";

MergeUMIsStrategySimple::MergeUMIsStrategySimple(unsigned int max_merge_distance)
		: _max_merge_distance(max_merge_distance)
{
	srand(42);
}

void MergeUMIsStrategySimple::merge(CellsDataContainer &container) const
{
	Tools::trace_time("UMI merge start");
	size_t total_cell_merged = 0, total_umi_merged = 0, total_cells_processed = 0;
	for (size_t cell_id = 0; cell_id < container.total_cells_number(); ++cell_id)
	{
		auto const &cell = container.cell(cell_id);
		if (!cell.is_real())
			continue;

		total_cells_processed++;
		for (auto const &gene : cell.genes())
		{
			s_hash_t bad_umis;
			for (auto const &umi : gene.second.umis())
			{
				if (this->is_umi_real(umi.first))
					continue;

				bad_umis.insert(umi.first);
			}

			if (bad_umis.empty())
				continue;

			CellsDataContainer::s_s_hash_t merge_targets = this->find_targets(gene.second.umis(), bad_umis);

			total_cell_merged++;
			total_umi_merged += merge_targets.size();

			container.merge_umis(cell_id, gene.first, merge_targets);
		}
	}

	L_TRACE << total_cells_processed << " cells processed. Merged " << total_umi_merged << " UMIs from "
	        << total_cell_merged << " cells.";
	Tools::trace_time("UMI merge finished");
}

bool MergeUMIsStrategySimple::is_umi_real(const std::string &umi) const
{
	return umi.find('N') == std::string::npos;
}

CellsDataContainer::s_s_hash_t MergeUMIsStrategySimple::find_targets(const Gene::umis_t &all_umis,
                                                                     const s_hash_t &bad_umis) const
{
	s_vec_t umis_without_pairs;
	CellsDataContainer::s_s_hash_t merge_targets;
	for (auto const &bad_umi : bad_umis)
	{
		int min_ed = std::numeric_limits<unsigned>::max();
		std::string best_target = "";
		long best_target_size = 0;
		for (auto const &target_umi : all_umis)
		{
			if (bad_umis.find(target_umi.first) != bad_umis.end())
				continue;

			unsigned ed = Tools::hamming_distance(target_umi.first, bad_umi);
			if (ed < min_ed || (ed == min_ed && target_umi.second.read_count > best_target_size))
			{
				min_ed = ed;
				best_target = target_umi.first;
				best_target_size = target_umi.second.read_count;
			}
		}

		if (best_target == "" || min_ed > this->_max_merge_distance)
		{
			umis_without_pairs.push_back(bad_umi);
		}
		else
		{
			merge_targets[bad_umi] = best_target;
		}
	}

	if (!umis_without_pairs.empty())
	{
		this->remove_similar_wrong_umis(umis_without_pairs);
		auto merge_targets2 = this->fill_wrong_umis(umis_without_pairs);
		merge_targets.insert(merge_targets2.begin(), merge_targets2.end());
	}

	return merge_targets;
}

void MergeUMIsStrategySimple::remove_similar_wrong_umis(s_vec_t &wrong_umis) const
{
	s_vec_t res;
	for (size_t i = 0; i < wrong_umis.size(); ++i)
	{
		for (size_t j = i + 1; j < wrong_umis.size(); ++j)
		{
			if (Tools::hamming_distance(wrong_umis[i], wrong_umis[j]) == 0)
				goto skip_umi;
		}

		res.push_back(wrong_umis[i]);
		skip_umi:;
	}
	wrong_umis = res;
}

CellsDataContainer::s_s_hash_t MergeUMIsStrategySimple::fill_wrong_umis(s_vec_t &wrong_umis) const
{
	CellsDataContainer::s_s_hash_t merge_targets;
	for (auto &umi : wrong_umis)
	{
		std::string target_umi(umi);
		for (std::string::size_type i = 0; i < target_umi.size(); ++i)
		{
			if (target_umi[i] != 'N')
				continue;

			target_umi[i] = MergeUMIsStrategySimple::nucleotides[rand() % MergeUMIsStrategySimple::nucleotides.size()];
		}
		merge_targets[umi] = target_umi;
	}
	return merge_targets;
}
}
}
}