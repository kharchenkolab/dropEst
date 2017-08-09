#include "MergeStrategyFactory.h"

#include "MergeStrategyAbstract.h"
#include "BrokenRealBarcodesMergeStrategy.h"
#include "DummyMergeStrategy.h"
#include "PoissonRealBarcodesMergeStrategy.h"
#include "RealBarcodesMergeStrategy.h"
#include "SimpleMergeStrategy.h"
#include "PoissonSimpleMergeStrategy.h"

#include <Estimation/Merge/BarcodesParsing/InDropBarcodesParser.h>
#include <Estimation/Merge/BarcodesParsing/ConstLengthBarcodesParser.h>
#include <Estimation/Merge/UMIs/MergeUMIsStrategySimple.h>

#include <fstream>

namespace Estimation
{
namespace Merge
{
	MergeStrategyFactory::MergeStrategyFactory(const boost::property_tree::ptree &config, int min_genes_after_merge)
	{
		auto main_config = config.get_child("Merge", boost::property_tree::ptree());

		this->_min_genes_before_merge = main_config.get<size_t>("min_genes_before_merge", 10);

		if (min_genes_after_merge > 0)
		{
			this->_min_genes_after_merge = unsigned(min_genes_after_merge);
		}
		else
		{
			this->_min_genes_after_merge = main_config.get<size_t>("min_genes_after_merge", 10);
		}

		this->_max_merge_edit_distance = main_config.get<unsigned>("max_cb_merge_edit_distance");

		this->_min_merge_fraction = main_config.get<double>("min_merge_fraction", 0.2);

		this->_barcodes_type = main_config.get<std::string>("barcodes_type", "indrop");
		this->_barcodes_filename = main_config.get<std::string>("barcodes_file", "");
		if (!this->_barcodes_filename.empty() && !std::ifstream(this->_barcodes_filename))
			throw std::runtime_error("Can't open barcodes file '" + this->_barcodes_filename + "'");

		auto poisson_config = config.get_child("PreciseMerge", boost::property_tree::ptree());

		this->_max_merge_prob = poisson_config .get<double>("max_merge_prob", 1e-4);
		this->_max_real_cb_merge_prob = poisson_config.get<double>("max_real_merge_prob", 1e-7);

		this->_max_umi_merge_edit_distance = main_config.get<unsigned>("max_umi_merge_edit_distance", 2);
	}

	MergeStrategyFactory::merge_cb_ptr MergeStrategyFactory::get_cb_strat(bool merge_tags, bool use_poisson) const
	{
		if (!merge_tags)
			return merge_cb_ptr(new DummyMergeStrategy(this->_min_genes_before_merge, this->_min_genes_after_merge));

		if (!use_poisson)
			return this->get_cb_strat();

		return this->get_cb_poisson_strat();
	}

	MergeStrategyFactory::merge_cb_ptr MergeStrategyFactory::get_cb_strat() const
	{
		if (this->_barcodes_filename.empty())
			return merge_cb_ptr(new SimpleMergeStrategy(this->_min_genes_before_merge, this->_min_genes_after_merge,
			                                            this->_max_merge_edit_distance, this->_min_merge_fraction));

		return merge_cb_ptr(new RealBarcodesMergeStrategy(this->get_barcodes_parser(),
		                                                  this->_min_genes_before_merge, this->_min_genes_after_merge,
		                                                  this->_max_merge_edit_distance, this->_min_merge_fraction));

//		if (merge_type == "broken")
//			return std::shared_ptr<MergeStrategyAbstract>(new BrokenRealBarcodesMergeStrategy(barcodes_filename, merge_config));
	}

	MergeStrategyFactory::merge_cb_ptr MergeStrategyFactory::get_cb_poisson_strat() const
	{
		PoissonTargetEstimator target_estimator(this->_max_merge_prob, this->_max_real_cb_merge_prob);

		if (this->_barcodes_filename.empty())
			return merge_cb_ptr(new PoissonSimpleMergeStrategy(target_estimator, this->_min_genes_before_merge,
			                                                   this->_min_genes_after_merge,
			                                                   this->_max_merge_edit_distance));

		return merge_cb_ptr(new PoissonRealBarcodesMergeStrategy(target_estimator, this->get_barcodes_parser(),
		                                                         this->_min_genes_before_merge,
		                                                         this->_min_genes_after_merge,
		                                                         this->_max_merge_edit_distance));
	}

	MergeStrategyFactory::merge_umi_ptr MergeStrategyFactory::get_umi() const
	{
		return merge_umi_ptr(new UMIs::MergeUMIsStrategySimple(this->_max_umi_merge_edit_distance));
	}

	MergeStrategyFactory::barcodes_parser_ptr MergeStrategyFactory::get_barcodes_parser() const
	{
		using namespace BarcodesParsing;
		if (this->_barcodes_type.empty())
			throw std::runtime_error("Empty barcodes type!");

		if (this->_barcodes_type == "indrop")
			return barcodes_parser_ptr(new InDropBarcodesParser(this->_barcodes_filename));

		if (this->_barcodes_type == "const")
			return barcodes_parser_ptr(new ConstLengthBarcodesParser(this->_barcodes_filename));

		throw std::runtime_error("Unexpected barcodes type: " + this->_barcodes_type);
	}
}
}