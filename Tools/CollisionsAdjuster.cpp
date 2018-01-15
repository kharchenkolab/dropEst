#include <cmath>
#include "CollisionsAdjuster.h"
#include "UtilFunctions.h"

namespace Tools
{
	CollisionsAdjuster::CollisionsAdjuster()
		: _sum_collisions(0)
		, _last_total_gene_size(0)
	{}

	void CollisionsAdjuster::init(const probs_vec_t &umi_probabilities, size_t max_gene_expression)
	{
		this->_sum_collisions = 0;
		this->_last_total_gene_size = 0;
		this->_umi_probabilities = umi_probabilities;
		this->_umi_probabilities_neg_prod = probs_vec_t(umi_probabilities.size(), 1);
		CollisionsAdjuster::update_adjusted_sizes(max_gene_expression);
	}

	void CollisionsAdjuster::update_adjusted_sizes(size_t max_gene_expression)
	{
		for (size_t s = this->_adjusted_sizes.size() + 1; s <= max_gene_expression; ++s)
		{
			const size_t total_size = s + size_t(this->_sum_collisions);
			double new_umi_prob = 0;
			for (int i = 0; i < this->_umi_probabilities.size(); ++i)
			{
				this->_umi_probabilities_neg_prod[i] *= fpow(1 - this->_umi_probabilities[i], total_size - this->_last_total_gene_size);
				new_umi_prob += this->_umi_probabilities[i] * (1 - this->_umi_probabilities_neg_prod[i]);
			}

			this->_last_total_gene_size = total_size;

			const double collision_num = 1.0 / (1.0 - new_umi_prob) - 1.0;
			this->_sum_collisions += collision_num;
			this->_adjusted_sizes.push_back(std::lround(s + this->_sum_collisions));
		}
	}

	size_t CollisionsAdjuster::estimate_adjusted_gene_expression(size_t expression)
	{
		if (expression > this->_adjusted_sizes.size())
		{
			this->update_adjusted_sizes(expression);
		}

		return this->_adjusted_sizes.at(expression - 1);
	}
}