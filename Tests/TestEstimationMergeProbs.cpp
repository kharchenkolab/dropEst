#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>
#include <boost/unordered_map.hpp>

#include "Estimation/CellsDataContainer.h"
#include "Tools/Logs.h"

#include <RInside.h>
#include <Estimation/Merge/MergeStrategyFactory.h>
#include <Estimation/Merge/SimpleMergeStrategy.h>
#include <Estimation/Merge/RealBarcodesMergeStrategy.h>
#include <Estimation/Merge/PoissonRealBarcodesMergeStrategy.h>

using namespace Estimation;

struct Fixture
{
	Fixture()
		: real_cb_strat(new Merge::PoissonRealBarcodesMergeStrategy(PROJ_DATA_PATH + std::string("/barcodes.txt"), 9, 0, 0, 100, 0))
		, container_full(real_cb_strat, 1)
	{
		Tools::init_test_logs(boost::log::trivial::info);
		CellsDataContainer::s_i_map_t cells_ids;
		this->container_full.add_record("AAATTAGGTCCA", "AAACCT", "Gene1", cells_ids);
		this->container_full.add_record("AAATTAGGTCCA", "CCCCCT", "Gene2", cells_ids);
		this->container_full.add_record("AAATTAGGTCCA", "ACCCCT", "Gene3", cells_ids);

		this->container_full.add_record("AAATTAGGTCCC", "CAACCT", "Gene1", cells_ids);

		this->container_full.add_record("AAATTAGGTCCG", "CAACCT", "Gene1", cells_ids);

		this->container_full.add_record("AAATTAGGTCGG", "AAACCT", "Gene1", cells_ids);
		this->container_full.add_record("AAATTAGGTCGG", "CCCCCT", "Gene2", cells_ids);

		this->container_full.add_record("CCCTTAGGTCCA", "CCATTC", "Gene3", cells_ids);
		this->container_full.add_record("CCCTTAGGTCCA", "CCCCCT", "Gene2", cells_ids);
		this->container_full.add_record("CCCTTAGGTCCA", "ACCCCT", "Gene3", cells_ids);

		this->container_full.add_record("CAATTAGGTCCG", "CAACCT", "Gene1", cells_ids);
		this->container_full.add_record("CAATTAGGTCCG", "AAACCT", "Gene1", cells_ids);
		this->container_full.add_record("CAATTAGGTCCG", "CCCCCT", "Gene2", cells_ids);
		this->container_full.add_record("CAATTAGGTCCG", "TTTTTT", "Gene2", cells_ids);
		this->container_full.add_record("CAATTAGGTCCG", "TTCTTT", "Gene2", cells_ids);

		this->container_full.add_record("CCCCCCCCCCCC", "CAACCT", "Gene1", cells_ids);
		this->container_full.add_record("CCCCCCCCCCCC", "AAACCT", "Gene1", cells_ids);
		this->container_full.add_record("CCCCCCCCCCCC", "CCCCCT", "Gene2", cells_ids);
		this->container_full.add_record("CCCCCCCCCCCC", "TTTTTT", "Gene2", cells_ids);
		this->container_full.add_record("CCCCCCCCCCCC", "TTCTTT", "Gene2", cells_ids);
		this->container_full.set_initialized();
	}

	std::shared_ptr<Merge::PoissonRealBarcodesMergeStrategy> real_cb_strat;

	CellsDataContainer container_full;
};

BOOST_AUTO_TEST_SUITE(TestEstimatorMergeProbs)

	BOOST_FIXTURE_TEST_CASE(testPoissonMergeInit, Fixture)
	{
		this->real_cb_strat->init(this->container_full);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution.size(), 7);

		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["AAACCT"], 4);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["CCCCCT"], 5);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["ACCCCT"], 2);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["CAACCT"], 4);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["CCATTC"], 1);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["TTTTTT"], 2);
		BOOST_CHECK_EQUAL(this->real_cb_strat->_umis_distribution["TTCTTT"], 2);

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(5).size(), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(6).size(), 2);

		std::map<std::string, int> distr;
		for (auto const &umi : this->real_cb_strat->_umis_bootstrap_distribution)
		{
			distr[umi]++;
		}

		for (auto const &umi : distr)
		{
			BOOST_CHECK_EQUAL(umi.second, this->real_cb_strat->_umis_distribution[umi.first]);
		}

	}

	BOOST_FIXTURE_TEST_CASE(testPoissonMergeProbs, Fixture)
	{
		this->real_cb_strat->init(this->container_full);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_bootstrap_intersect_prob(this->container_full, 0, 1), 1);
		BOOST_CHECK_LE(std::abs(this->real_cb_strat->get_bootstrap_intersect_prob(this->container_full, 1, 2) - 0.16), 0.05);
		BOOST_CHECK_LE(std::abs(this->real_cb_strat->get_bootstrap_intersect_prob(this->container_full, 3, 4) - 0.15), 0.05);
		BOOST_CHECK_LE(std::abs(this->real_cb_strat->get_bootstrap_intersect_prob(this->container_full, 5, 6, 100000) - 0.045), 0.01);
	}

BOOST_AUTO_TEST_SUITE_END()