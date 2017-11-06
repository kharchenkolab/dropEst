#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Estimation/CellsDataContainer.h"
#include "Tools/Logs.h"

#include <RInside.h>
#include <Estimation/Merge/MergeStrategyFactory.h>
#include <Estimation/Merge/SimpleMergeStrategy.h>
#include <Estimation/Merge/RealBarcodesMergeStrategy.h>
#include <Estimation/Merge/PoissonRealBarcodesMergeStrategy.h>
#include <Estimation/Merge/PoissonTargetEstimator.h>
#include <Estimation/Merge/UMIs/MergeUMIsStrategySimple.h>
#include <Estimation/Merge/BarcodesParsing/InDropBarcodesParser.h>

using namespace Estimation;

struct Fixture
{
	Fixture()
		: estimator(1.0e-4, 1.0e-7)
	{
		std::stringstream config("<Estimation>\n"
								 "        <barcodes_type>indrop</barcodes_type>\n"
								 "        <min_merge_fraction>0</min_merge_fraction>\n"
								 "        <max_merge_edit_distance>100</max_merge_edit_distance>\n"
								 "        <min_genes_after_merge>0</min_genes_after_merge>\n"
								 "        <min_genes_before_merge>0</min_genes_before_merge>\n"
								 "    </Estimation>");

		using Mark = UMI::Mark;

		auto barcodes_parser = std::shared_ptr<Merge::BarcodesParsing::BarcodesParser>(
				new Merge::BarcodesParsing::InDropBarcodesParser(PROJ_DATA_PATH + std::string("/barcodes/test_est")));

		const Merge::PoissonTargetEstimator target_estimator(1e-4, 1e-7);

		this->real_cb_strat = std::make_shared<Merge::PoissonRealBarcodesMergeStrategy>(target_estimator, barcodes_parser, 0, 0, 7);

		this->container_full = std::make_shared<CellsDataContainer>(this->real_cb_strat, std::make_shared<Merge::UMIs::MergeUMIsStrategySimple>(1),
		                                                            Mark::get_by_code(Mark::DEFAULT_CODE), -1);

		Tools::init_test_logs(boost::log::trivial::info);
		this->container_full->add_record("AAATTAGGTCCA", "AAACCT", "Gene1");
		this->container_full->add_record("AAATTAGGTCCA", "CCCCCT", "Gene2");
		this->container_full->add_record("AAATTAGGTCCA", "ACCCCT", "Gene3");

		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene1");

		this->container_full->add_record("AAATTAGGTCCG", "CAACCT", "Gene1");

		this->container_full->add_record("AAATTAGGTCGG", "AAACCT", "Gene1");
		this->container_full->add_record("AAATTAGGTCGG", "CCCCCT", "Gene2");

		this->container_full->add_record("CCCTTAGGTCCA", "CCATTC", "Gene3");
		this->container_full->add_record("CCCTTAGGTCCA", "CCCCCT", "Gene2");
		this->container_full->add_record("CCCTTAGGTCCA", "ACCCCT", "Gene3");

		this->container_full->add_record("CAATTAGGTCCG", "CAACCT", "Gene1");
		this->container_full->add_record("CAATTAGGTCCG", "AAACCT", "Gene1");
		this->container_full->add_record("CAATTAGGTCCG", "CCCCCT", "Gene2");
		this->container_full->add_record("CAATTAGGTCCG", "TTTTTT", "Gene2");
		this->container_full->add_record("CAATTAGGTCCG", "TTCTTT", "Gene2");

		this->container_full->add_record("CCCCCCCCCCCC", "CAACCT", "Gene1");
		this->container_full->add_record("CCCCCCCCCCCC", "AAACCT", "Gene1");
		this->container_full->add_record("CCCCCCCCCCCC", "CCCCCT", "Gene2");
		this->container_full->add_record("CCCCCCCCCCCC", "TTTTTT", "Gene2");
		this->container_full->add_record("CCCCCCCCCCCC", "TTCTTT", "Gene2");

		this->container_full->add_record("TAATTAGGTCCA", "AAAAAA", "Gene4");
		this->container_full->set_initialized();
	}

	std::shared_ptr<Merge::PoissonRealBarcodesMergeStrategy> real_cb_strat;
	Merge::PoissonTargetEstimator estimator;

	std::shared_ptr<CellsDataContainer> container_full;
};

BOOST_AUTO_TEST_SUITE(TestEstimatorMergeProbs)

	BOOST_FIXTURE_TEST_CASE(testPoissonMergeInit, Fixture)
	{
		this->estimator.init(this->container_full->umi_distribution());
		BOOST_CHECK_EQUAL(this->estimator._umi_distribution.size(), 8);

		BOOST_CHECK_EQUAL(this->container_full->cell(5).size(), 2);
		BOOST_CHECK_EQUAL(this->container_full->cell(6).size(), 2);

//		std::map<Merge::PoissonRealBarcodesMergeStrategy::bs_umi_t, int> distr;
//		for (auto const &umi : this->real_cb_strat->_umis_bootstrap_distribution)
//		{
//			distr[umi]++;
//		}
//
//		for (auto const &umi : distr)
//		{
//			BOOST_CHECK_EQUAL(umi.second, this->real_cb_strat->_umis_distribution[umi.first]);
//		}
	}

	BOOST_FIXTURE_TEST_CASE(testIntersectionSizeEstimation, Fixture)
	{
		this->estimator.init(this->container_full->umi_distribution());

		// Values were obtained with R
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_genes_intersection_size(1, 5) - 0.7264), 1e-2);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_genes_intersection_size(2, 5) - 1.4484), 1e-2);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_genes_intersection_size(3, 5) - 2.1380), 1e-2);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_genes_intersection_size(4, 5) - 2.7923), 1e-2);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_genes_intersection_size(5, 5) - 3.4346), 1e-2);

		BOOST_CHECK_LE(std::abs(this->estimator.estimate_genes_intersection_size(5, 3) - 2.1380), 1e-2);
	}

	BOOST_FIXTURE_TEST_CASE(testPoissonMergeProbs, Fixture)
	{
		this->estimator.init(this->container_full->umi_distribution());
		BOOST_CHECK_EQUAL(this->estimator.estimate_intersection_prob(*this->container_full, 0, 1).merge_probability, 1);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_intersection_prob(*this->container_full, 1, 2).merge_probability - 0.16), 0.05);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_intersection_prob(*this->container_full, 3, 4).merge_probability - 0.15), 0.05);
		BOOST_CHECK_LE(std::abs(this->estimator.estimate_intersection_prob(*this->container_full, 5, 6).merge_probability - 0.05), 0.01);
	}

	BOOST_FIXTURE_TEST_CASE(testPoissonMergeRejections, Fixture)
	{
		this->real_cb_strat->init(*this->container_full);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 7), -1);
	}

BOOST_AUTO_TEST_SUITE_END()
