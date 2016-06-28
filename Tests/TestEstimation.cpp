#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>
#include <boost/unordered_map.hpp>

#include <sstream>
#include "Estimation/CellsDataContainer.h"
#include "Tools/Logs.h"

#ifdef R_LIBS
#include <RInside.h>
#endif

using namespace Estimation;

struct Fixture
{
	Fixture()
		: container(0, 0, 0, 0, 100, 10)
		, container_full(0, 0, true, 0, 100, 1)
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
		this->container_full.set_initialized();
	}

	CellsDataContainer container;
	CellsDataContainer container_full;
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
		CellsDataContainer::ISIHM cb_reassigned_to;
		CellsDataContainer::ints_t cb_reassigned;
		cb_reassigned.push_back(0);
		cb_reassigned.push_back(1);
		cb_reassigned.push_back(2);

		std::string cb1 = "CCCCCCCCCCCCC";
		std::string cb2 = "CCCCCCCCCCCCA";
		std::string cb3 = "CCCCCCCCCCCAA";
		container._cells_barcodes.push_back(cb1);
		container._cells_barcodes.push_back(cb2);
		container._cells_barcodes.push_back(cb3);

		container._cells_genes.resize(2);
		container._cells_genes[0]["G1"]["UMI1"] = 1;
		container._cells_genes[0]["G1"]["UMI2"] = 1;
		container._cells_genes[0]["G2"]["UMI1"] = 1;

		container._cells_genes[1]["G1"]["UMI1"] = 3;
		container._cells_genes[1]["G1"]["UMI3"] = 5;
		container._cells_genes[1]["G2"]["UMI2"] = 7;

		container.merge(1, 1, CellsDataContainer::IndexedValue(0, 3), cb_reassigned, cb_reassigned_to); //merge 0 to 1

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

	BOOST_FIXTURE_TEST_CASE(testBarcodesFile, Fixture)
	{
		std::vector<std::string> cb1;
		std::vector<std::string> cb2;

		CellsDataContainer::get_barcodes_list(PROJ_DATA_PATH + std::string("/barcodes.txt"), cb1, cb2);

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

		BOOST_CHECK_THROW(CellsDataContainer::get_barcodes_list("/barcodes.txt", cb1, cb2), std::runtime_error);
	}

	BOOST_FIXTURE_TEST_CASE(testIncrement, Fixture)
	{
		boost::unordered_map<std::string, int> map;
		map["chr1"]++;
		BOOST_CHECK_EQUAL(map["chr1"], 1);
	}

	BOOST_FIXTURE_TEST_CASE(testUmigsIntersection, Fixture)
	{
		size_t is1 = this->container_full.get_umigs_intersection_size(this->container_full._cell_ids_by_cb["AAATTAGGTCCA"],
		                                                              this->container_full._cell_ids_by_cb["CCCTTAGGTCCA"]);

		size_t is2 = this->container_full.get_umigs_intersection_size(this->container_full._cell_ids_by_cb["AAATTAGGTCCC"],
		                                                              this->container_full._cell_ids_by_cb["AAATTAGGTCCG"]);

		size_t is3 = this->container_full.get_umigs_intersection_size(this->container_full._cell_ids_by_cb["AAATTAGGTCCA"],
		                                                              this->container_full._cell_ids_by_cb["AAATTAGGTCCC"]);

		BOOST_CHECK_EQUAL(is1, 2);
		BOOST_CHECK_EQUAL(is2, 1);
		BOOST_CHECK_EQUAL(is3, 0);
	}

	BOOST_FIXTURE_TEST_CASE(testFillDistances, Fixture)
	{
		static const std::string ar_cbs1[] = {"AAT", "AAA", "CCT"};
		std::vector<std::string> cbs1(ar_cbs1, ar_cbs1 + sizeof(ar_cbs1) / sizeof(ar_cbs1[0]));

		CellsDataContainer::i_counter_t d1, d2;
		CellsDataContainer::fill_distances_to_cb("ACT", "ACT", cbs1, cbs1, d1, d2);

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
		CellsDataContainer::fill_distances_to_cb("CAA", "TTAGGTCCG", cbs1, cbs2, d1, d2);

		CellsDataContainer::ids_t ids = this->container_full.get_real_neighbour_cbs(cbs1, cbs2, "CAATTAGGTCCG", d1, d2);

		BOOST_CHECK_EQUAL(ids.size(), 2);

		if (ids.size() == 2)
		{
			BOOST_CHECK_EQUAL(this->container_full._cells_barcodes[ids[0]], "AAATTAGGTCCA");
			BOOST_CHECK_EQUAL(this->container_full._cells_barcodes[ids[1]], "AAATTAGGTCCC");
		}

		d1.clear(); d2.clear();
		CellsDataContainer::fill_distances_to_cb("AAA", "TTAGGTCCC", cbs1, cbs2, d1, d2);
		ids = this->container_full.get_real_neighbour_cbs(cbs1, cbs2, "CAATTAGGTCCG", d1, d2);

		BOOST_CHECK_EQUAL(ids.size(), 1);

		if (ids.size() >= 2)
		{
			BOOST_CHECK_EQUAL(this->container_full._cells_barcodes[ids[0]], "AAATTAGGTCCC");
		}
	}

	BOOST_FIXTURE_TEST_CASE(testRealNeighbours, Fixture)
	{
		static const std::string ar_cbs1[] = {"AAT", "GAA", "AAA"};
		static const std::string ar_cbs2[] = {"TTAGGTCCA", "TTAGGGGCC", "TTAGGTCCC"};
		CellsDataContainer::names_t cbs1(ar_cbs1, ar_cbs1 + sizeof(ar_cbs1) / sizeof(ar_cbs1[0]));
		CellsDataContainer::names_t cbs2(ar_cbs2, ar_cbs2 + sizeof(ar_cbs2) / sizeof(ar_cbs2[0]));

		BOOST_CHECK_EQUAL(this->container_full.get_real_cb(0, cbs1, cbs2, 9), 0);
		BOOST_CHECK_EQUAL(this->container_full.get_real_cb(1, cbs1, cbs2, 9), 1);
		BOOST_CHECK_EQUAL(this->container_full.get_real_cb(2, cbs1, cbs2, 9), 1);
		BOOST_CHECK_EQUAL(this->container_full.get_real_cb(3, cbs1, cbs2, 9), 0);
		BOOST_CHECK_EQUAL(this->container_full.get_real_cb(4, cbs1, cbs2, 9), 0);
		BOOST_CHECK_EQUAL(this->container_full.get_real_cb(5, cbs1, cbs2, 9), 0);
	}

	BOOST_FIXTURE_TEST_CASE(testMergeByRealBarcodes, Fixture)
	{
		this->container_full.merge_by_real_barcodes(PROJ_DATA_PATH + std::string("/barcodes.txt"), 9);

		BOOST_CHECK_EQUAL(this->container_full.filtered_cells().size(), 2);

		if (this->container_full.filtered_cells().size() < 2)
			return;

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).size(), 3);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).size(), 1);

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene1").size(), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene1").at("AAACCT"), 3);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene2").size(), 1);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene2").at("CCCCCT"), 4);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene3").size(), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene3").at("ACCCCT"), 2);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene3").at("CCATTC"), 1);

		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[1]).at("Gene1").size(), 1);
		BOOST_CHECK_EQUAL(this->container_full.cell_genes(this->container_full.filtered_cells()[0]).at("Gene1").at("CAACCT"), 2);
	}

BOOST_AUTO_TEST_SUITE_END()