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

using namespace Estimation;

struct Fixture
{
	Fixture()
		: base_strat(new Merge::SimpleMergeStrategy(0, 0, 100, 0))
		, real_cb_strat(new Merge::RealBarcodesMergeStrategy(PROJ_DATA_PATH + std::string("/barcodes.txt"), 9, 0, 0, 7, 0))
		, container(base_strat, 10)
		, container_full(real_cb_strat, 1)
	{
		Tools::init_test_logs(boost::log::trivial::info); //TODO Fix tests!!!
		CellsDataContainer::s_i_map_t cells_ids;
		this->container_full.add_record("AAATTAGGTCCA", "AAACCT", "Gene1", cells_ids); //0, real
		this->container_full.add_record("AAATTAGGTCCA", "CCCCCT", "Gene2", cells_ids);
		this->container_full.add_record("AAATTAGGTCCA", "ACCCCT", "Gene3", cells_ids);
		this->container_full.add_record("AAATTAGGTCCA", "ACCCCT", "Gene4", cells_ids);

		this->container_full.add_record("AAATTAGGTCCC", "CAACCT", "Gene1", cells_ids); //1, real
		this->container_full.add_record("AAATTAGGTCCC", "CAACCT", "Gene10", cells_ids);
		this->container_full.add_record("AAATTAGGTCCC", "CAACCT", "Gene20", cells_ids);

		this->container_full.add_record("AAATTAGGTCCG", "CAACCT", "Gene1", cells_ids); //2, false

		this->container_full.add_record("AAATTAGGTCGG", "AAACCT", "Gene1", cells_ids); //3, false
		this->container_full.add_record("AAATTAGGTCGG", "CCCCCT", "Gene2", cells_ids);

		this->container_full.add_record("CCCTTAGGTCCA", "CCATTC", "Gene3", cells_ids); //4, false
		this->container_full.add_record("CCCTTAGGTCCA", "CCCCCT", "Gene2", cells_ids);
		this->container_full.add_record("CCCTTAGGTCCA", "ACCCCT", "Gene3", cells_ids);

		this->container_full.add_record("CAATTAGGTCCG", "CAACCT", "Gene1", cells_ids); //5, false
		this->container_full.add_record("CAATTAGGTCCG", "AAACCT", "Gene1", cells_ids);
		this->container_full.add_record("CAATTAGGTCCG", "CCCCCT", "Gene2", cells_ids);

		this->container_full.add_record("AAAAAAAAAAAA", "CCCCCT", "Gene2", cells_ids); //6, false, excluded
		this->container_full.set_initialized();
	}

	std::shared_ptr<Merge::SimpleMergeStrategy> base_strat;
	std::shared_ptr<Merge::RealBarcodesMergeStrategy> real_cb_strat;

	CellsDataContainer container;
	CellsDataContainer container_full;
};

BOOST_AUTO_TEST_SUITE(TestEstimator)

	BOOST_FIXTURE_TEST_CASE(testR, Fixture)
	{

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
	}

	BOOST_FIXTURE_TEST_CASE(testBarcodesFile, Fixture)
	{
		std::vector<std::string> cb1;
		std::vector<std::string> cb2;

		Merge::RealBarcodesMergeStrategy::get_barcodes_list(PROJ_DATA_PATH + std::string("/barcodes.txt"), cb1, cb2);

		BOOST_CHECK_EQUAL(cb1.size(), 3);
		BOOST_CHECK_EQUAL(cb2.size(), 3);

		if (cb1.size() >= 3)
		{
			BOOST_CHECK_EQUAL(cb1[0], "AAT");
			BOOST_CHECK_EQUAL(cb1[1], "GAA");
			BOOST_CHECK_EQUAL(cb1[2], "AAA");
		}

		if (cb2.size() >= 3)
		{
			BOOST_CHECK_EQUAL(cb2[0], "TTAGGTCCA");
			BOOST_CHECK_EQUAL(cb2[1], "TTAGGGGCC");
			BOOST_CHECK_EQUAL(cb2[2], "TTAGGTCCC");
		}

		BOOST_CHECK_THROW(Merge::RealBarcodesMergeStrategy::get_barcodes_list("/barcodes.txt", cb1, cb2), std::runtime_error);
	}

	BOOST_FIXTURE_TEST_CASE(testIncrement, Fixture)
	{
		boost::unordered_map<std::string, int> map;
		map["chr1"]++;
		BOOST_CHECK_EQUAL(map["chr1"], 1);
		map.emplace(std::make_pair("chr1", 0));
		BOOST_CHECK_EQUAL(map["chr1"], 1);
	}

	BOOST_FIXTURE_TEST_CASE(testUmigsIntersection, Fixture)
	{
		using namespace Merge;
		auto cell_ids_by_cb = this->container_full.cell_ids_by_cb();
		long is1 = RealBarcodesMergeStrategy::get_umigs_intersect_size(this->container_full.cell_genes(cell_ids_by_cb["AAATTAGGTCCA"]),
																		 this->container_full.cell_genes(cell_ids_by_cb["CCCTTAGGTCCA"]));

		long is2 = RealBarcodesMergeStrategy::get_umigs_intersect_size(this->container_full.cell_genes(cell_ids_by_cb["AAATTAGGTCCC"]),
																		 this->container_full.cell_genes(cell_ids_by_cb["AAATTAGGTCCG"]));

		long is3 = RealBarcodesMergeStrategy::get_umigs_intersect_size(this->container_full.cell_genes(cell_ids_by_cb["AAATTAGGTCCA"]),
																		 this->container_full.cell_genes(cell_ids_by_cb["AAATTAGGTCCC"]));

		BOOST_CHECK_EQUAL(is1, 2);
		BOOST_CHECK_EQUAL(is2, 1);
		BOOST_CHECK_EQUAL(is3, 0);
	}

	BOOST_FIXTURE_TEST_CASE(testFillDistances, Fixture)
	{
		static const std::string ar_cbs1[] = {"AAT", "AAA", "CCT"};
		std::vector<std::string> cbs1(ar_cbs1, ar_cbs1 + sizeof(ar_cbs1) / sizeof(ar_cbs1[0]));

		CellsDataContainer::i_counter_t d1, d2;

		this->real_cb_strat->_barcodes1 = cbs1;
		this->real_cb_strat->_barcodes2 = cbs1;
		this->real_cb_strat->fill_distances_to_cb("ACT", "ACT", d1, d2);

		BOOST_CHECK_EQUAL(d1.size(), 3);
		BOOST_CHECK_EQUAL(d2.size(), 3);

		BOOST_CHECK_EQUAL(d1[0].value, 1);
		BOOST_CHECK_EQUAL(d1[1].value, 1);
		BOOST_CHECK_EQUAL(d1[2].value, 2);
		BOOST_CHECK_EQUAL(d1[2].index, 1);

		BOOST_CHECK_EQUAL(d2[0].value, 1);
		BOOST_CHECK_EQUAL(d2[1].value, 1);
		BOOST_CHECK_EQUAL(d2[2].value, 2);
		BOOST_CHECK_EQUAL(d2[2].index, 1);
	}

	BOOST_FIXTURE_TEST_CASE(testRealNeighboursCbs, Fixture)
	{
		static const std::string ar_cbs1[] = {"AAT", "GAA", "AAA"};
		static const std::string ar_cbs2[] = {"TTAGGTCCA", "TTAGGGGCC", "TTAGGTCCC"};
		CellsDataContainer::names_t cbs1(ar_cbs1, ar_cbs1 + sizeof(ar_cbs1) / sizeof(ar_cbs1[0]));
		CellsDataContainer::names_t cbs2(ar_cbs2, ar_cbs2 + sizeof(ar_cbs2) / sizeof(ar_cbs2[0]));

		CellsDataContainer::i_counter_t d1, d2;
		this->real_cb_strat->_barcodes1 = cbs1;
		this->real_cb_strat->_barcodes2 = cbs2;
		this->real_cb_strat->fill_distances_to_cb("CAA", "TTAGGTCCG", d1, d2);

		CellsDataContainer::ids_t ids = this->real_cb_strat->get_real_neighbour_cbs(this->container_full, this->container_full.cell_ids_by_cb().at("CAATTAGGTCCG"), d1, d2);

		BOOST_REQUIRE_EQUAL(ids.size(), 2);

		BOOST_CHECK_EQUAL(this->container_full.cell_barcode(ids[0]), "AAATTAGGTCCA");
		BOOST_CHECK_EQUAL(this->container_full.cell_barcode(ids[1]), "AAATTAGGTCCC");

		d1.clear(); d2.clear();
		this->real_cb_strat->fill_distances_to_cb("AAA", "TTAGGTCCC", d1, d2);
		ids = this->real_cb_strat->get_real_neighbour_cbs(this->container_full, this->container_full.cell_ids_by_cb().at("CAATTAGGTCCG"), d1, d2);

		BOOST_CHECK_EQUAL(ids.size(), 1);

		if (ids.size() >= 2)
		{
			BOOST_CHECK_EQUAL(this->container_full.cell_barcode(ids[0]), "AAATTAGGTCCC");
		}
	}

	BOOST_FIXTURE_TEST_CASE(testRealNeighbours, Fixture)
	{
		static const std::string ar_cbs1[] = {"AAT", "GAA", "AAA"};
		static const std::string ar_cbs2[] = {"TTAGGTCCA", "TTAGGGGCC", "TTAGGTCCC"};
		CellsDataContainer::names_t cbs1(ar_cbs1, ar_cbs1 + sizeof(ar_cbs1) / sizeof(ar_cbs1[0]));
		CellsDataContainer::names_t cbs2(ar_cbs2, ar_cbs2 + sizeof(ar_cbs2) / sizeof(ar_cbs2[0]));
		this->real_cb_strat->_barcodes1 = cbs1;
		this->real_cb_strat->_barcodes2 = cbs2;

		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(this->container_full, 0), 0);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(this->container_full, 1), 1);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(this->container_full, 2), 1);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(this->container_full, 3), 0);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(this->container_full, 4), 0);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(this->container_full, 5), 0);
	}

	BOOST_FIXTURE_TEST_CASE(testMergeByRealBarcodes, Fixture)
	{
		this->container_full.merge_and_filter();

		BOOST_CHECK_EQUAL(this->container_full.cell_barcodes_raw().size(), 7);
		BOOST_CHECK_EQUAL(this->container_full.filtered_cells().size(), 2);

		BOOST_REQUIRE(this->container_full.filtered_cells().size() >= 2);

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).size(), 3);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).size(), 4);

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene1").size(), 1);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene1").at("CAACCT"), 2);

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene1").size(), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene1").at("AAACCT"), 3);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene2").size(), 1);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene2").at("CCCCCT"), 4);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene3").size(), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene3").at("ACCCCT"), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene3").at("CCATTC"), 1);

		BOOST_CHECK(!this->container_full.is_cell_merged(0));
		BOOST_CHECK(!this->container_full.is_cell_merged(1));
		BOOST_CHECK(this->container_full.is_cell_merged(2));
		BOOST_CHECK(this->container_full.is_cell_merged(3));
		BOOST_CHECK(this->container_full.is_cell_merged(4));
		BOOST_CHECK(this->container_full.is_cell_merged(5));
		BOOST_CHECK(!this->container_full.is_cell_merged(6));

		BOOST_CHECK_EQUAL(this->container_full.excluded_cells().size(), 1);
	}

BOOST_AUTO_TEST_SUITE_END()