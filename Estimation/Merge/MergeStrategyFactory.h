#pragma once

#include <memory>
#include <string>

#include <boost/property_tree/ptree.hpp>

namespace Estimation
{
namespace Merge
{
	namespace BarcodesParsing
	{
		class BarcodesParser;
	}

	namespace UMIs
	{
		class MergeUMIsStrategySimple;
	}

	class MergeStrategyAbstract;

	class MergeStrategyFactory
	{
	public:
		typedef std::shared_ptr<MergeStrategyAbstract> merge_cb_ptr;
		typedef std::shared_ptr<UMIs::MergeUMIsStrategySimple> merge_umi_ptr;
		typedef std::shared_ptr<BarcodesParsing::BarcodesParser> barcodes_parser_ptr;

	private:
		std::string _merge_type;
		size_t _min_genes_before_merge;
		size_t _min_genes_after_merge;
		unsigned int _max_merge_edit_distance;

		double _min_merge_fraction;
		std::string _barcodes_filename;
		std::string _barcodes_type;

		double _max_merge_prob;
		double _max_real_cb_merge_prob;

		unsigned int _max_umi_merge_edit_distance;

	private:
		merge_cb_ptr get_cb_poisson_strat() const;
		merge_cb_ptr get_cb_strat() const;


		barcodes_parser_ptr get_barcodes_parser() const;

	public:
		merge_cb_ptr get_cb_strat(bool merge_tags, bool use_poisson) const;
		merge_umi_ptr get_umi() const;

		MergeStrategyFactory(const boost::property_tree::ptree &config, int min_genes_after_merge = -1);
	};
}
}