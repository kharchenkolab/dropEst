#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>

#include "Estimation/IndropResults.h"

#ifdef R_LIBS
#include <RInside.h>
#endif


#include <sstream>
#include <boost/unordered_map.hpp>

struct Fixture
{
	Fixture()
		: result()
	{}

	IndropResult result;
};

BOOST_AUTO_TEST_SUITE(TestTagsSearch)

	BOOST_FIXTURE_TEST_CASE(test1, Fixture)
	{
#ifdef R_LIBS
		using namespace Rcpp;
		boost::unordered_map<std::string, int> map;
		map["first"] = 1;
		map["second"] = 2;

		SEXP res = wrap(map.begin(), map.end());

		RInside R(0, 0);
		R["saved_vec"] = res;
		R.parseEvalQ("saveRDS(saved_vec, 'test_rds.rds')");
#endif
	}

BOOST_AUTO_TEST_SUITE_END()