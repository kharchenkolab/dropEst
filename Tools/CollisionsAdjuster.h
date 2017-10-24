#pragma once

#include <vector>
#include <cstddef>

namespace Tools
{
	class CollisionsAdjuster
	{
	public:
		typedef std::vector<size_t> size_vec_t;
		typedef std::vector<double> probs_vec_t;

	private:
		size_vec_t _adjusted_sizes;
		probs_vec_t _umi_probabilities;
		probs_vec_t _umi_probabilities_neg_prod;

		double _sum_collisions;
		size_t _last_total_gene_size;

	private:
		void update_adjusted_sizes(size_t max_gene_expression);

	public:
		CollisionsAdjuster();
		void init(const probs_vec_t &umi_probabilities, size_t max_gene_expression = 0);
		size_t estimate_adjusted_gene_expression(size_t expression);
	};
}

