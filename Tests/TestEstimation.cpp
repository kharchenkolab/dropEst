#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>

#include "Estimation/GenesContainer.h"

#ifdef R_LIBS
#include <RInside.h>
#endif


#include <sstream>
#include <boost/unordered_map.hpp>

struct Fixture
{
	Fixture()
		: container(0, 0, 0, 0, 100, 10)
	{}

	GenesContainer container;
};

BOOST_AUTO_TEST_SUITE(TestEstimator)

	BOOST_FIXTURE_TEST_CASE(testR, Fixture)
	{
#ifdef R_LIBS
//		using namespace Rcpp;
//		boost::unordered_map<std::string, int> map;
//		map["first"] = 1;
//		map["second"] = 2;
//
//		SEXP res = wrap(map.begin(), map.end());
//
//		RInside R(0, 0);
//		R["saved_vec"] = res;
//		R.parseEvalQ("saveRDS(saved_vec, 'test_rds.rds')");
#endif
	}

	BOOST_FIXTURE_TEST_CASE(testMerge, Fixture)
	{
		GenesContainer::ISIHM cb_reassigned_to;
		GenesContainer::ints_t cb_reassigned;
		cb_reassigned.push_back(0);
		cb_reassigned.push_back(1);
		cb_reassigned.push_back(2);

		std::string cb1 = "CCCCCCCCCCCCC";
		std::string cb2 = "CCCCCCCCCCCCA";
		std::string cb3 = "CCCCCCCCCCCAA";
		container._cells_names.push_back(cb1);
		container._cells_names.push_back(cb2);
		container._cells_names.push_back(cb3);

		container._cells_genes.resize(2);
		container._cells_genes[0]["G1"]["UMI1"] = 1;
		container._cells_genes[0]["G1"]["UMI2"] = 1;
		container._cells_genes[0]["G2"]["UMI1"] = 1;

		container._cells_genes[1]["G1"]["UMI1"] = 3;
		container._cells_genes[1]["G1"]["UMI3"] = 5;
		container._cells_genes[1]["G2"]["UMI2"] = 7;

		container.merge(1, 1, GenesContainer::IndexedCount(0, 3), cb_reassigned, cb_reassigned_to); //merge 0 to 1

		BOOST_CHECK_EQUAL(cb_reassigned[0], 1);
		BOOST_CHECK_EQUAL(cb_reassigned[1], 1);
		BOOST_CHECK_EQUAL(cb_reassigned[2], 2);

		BOOST_CHECK_EQUAL(cb_reassigned_to.size(), 1);
		BOOST_CHECK_EQUAL(cb_reassigned_to[1].size(), 1);
		BOOST_CHECK_EQUAL(cb_reassigned_to[1].count(0), 1);

		BOOST_CHECK_EQUAL(container._cells_genes[1]["G1"]["UMI1"], 4);
		BOOST_CHECK_EQUAL(container._cells_genes[1]["G1"]["UMI2"], 1);
		BOOST_CHECK_EQUAL(container._cells_genes[1]["G2"]["UMI1"], 1);
		BOOST_CHECK_EQUAL(container._cells_genes[1]["G1"]["UMI3"], 5);
		BOOST_CHECK_EQUAL(container._cells_genes[1]["G2"]["UMI2"], 7);
	}

	BOOST_FIXTURE_TEST_CASE(testIncrement, Fixture)
	{
		boost::unordered_map<std::string, int> map;
		map["chr1"]++;
		BOOST_CHECK_EQUAL(map["chr1"], 1);
	}

BOOST_AUTO_TEST_SUITE_END()