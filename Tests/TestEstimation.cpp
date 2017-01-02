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
#include <Estimation/Merge/BarcodesParsing/InDropBarcodesParser.h>
#include <Estimation/Merge/BarcodesParsing/ConstLengthBarcodesParser.h>

using namespace Estimation;

struct Fixture
{
	Fixture()
	{
		std::stringstream config("<Estimation>\n"
				"        <barcodes_type>indrop</barcodes_type>\n"
				"        <min_merge_fraction>0</min_merge_fraction>\n"
				"        <max_merge_edit_distance>7</max_merge_edit_distance>\n"
				"        <min_genes_after_merge>0</min_genes_after_merge>\n"
				"        <min_genes_before_merge>0</min_genes_before_merge>\n"
				"    </Estimation>");

		boost::property_tree::ptree pt;
		read_xml(config, pt);
		this->real_cb_strat = std::make_shared<Merge::RealBarcodesMergeStrategy>(PROJ_DATA_PATH + std::string("/barcodes/test_est"), pt.get_child("Estimation"));
		this->container_full = std::make_shared<CellsDataContainer>(this->real_cb_strat, 1);

		Tools::init_test_logs(boost::log::trivial::info);
		CellsDataContainer::s_i_map_t cells_ids;
		this->container_full->add_record("AAATTAGGTCCA", "AAACCT", "Gene1", cells_ids); //0, real
		this->container_full->add_record("AAATTAGGTCCA", "CCCCCT", "Gene2", cells_ids);
		this->container_full->add_record("AAATTAGGTCCA", "ACCCCT", "Gene3", cells_ids);
		this->container_full->add_record("AAATTAGGTCCA", "ACCCCT", "Gene4", cells_ids);

		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene1", cells_ids); //1, real
		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene10", cells_ids);
		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene20", cells_ids);

		this->container_full->add_record("AAATTAGGTCCG", "CAACCT", "Gene1", cells_ids); //2, false

		this->container_full->add_record("AAATTAGGTCGG", "AAACCT", "Gene1", cells_ids); //3, false
		this->container_full->add_record("AAATTAGGTCGG", "CCCCCT", "Gene2", cells_ids);

		this->container_full->add_record("CCCTTAGGTCCA", "CCATTC", "Gene3", cells_ids); //4, false
		this->container_full->add_record("CCCTTAGGTCCA", "CCCCCT", "Gene2", cells_ids);
		this->container_full->add_record("CCCTTAGGTCCA", "ACCCCT", "Gene3", cells_ids);

		this->container_full->add_record("CAATTAGGTCCG", "CAACCT", "Gene1", cells_ids); //5, false
		this->container_full->add_record("CAATTAGGTCCG", "AAACCT", "Gene1", cells_ids);
		this->container_full->add_record("CAATTAGGTCCG", "CCCCCT", "Gene2", cells_ids);

		this->container_full->add_record("AAAAAAAAAAAA", "CCCCCT", "Gene2", cells_ids); //6, false, excluded
		this->container_full->set_initialized();
	}

	std::shared_ptr<Merge::RealBarcodesMergeStrategy> real_cb_strat;
	std::shared_ptr<CellsDataContainer> container_full;
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
		auto cbs = this->real_cb_strat->_barcodes_parser->_barcodes;

		BOOST_CHECK_EQUAL(cbs[0].size(), 3);
		BOOST_CHECK_EQUAL(cbs[1].size(), 3);

		if (cbs[0].size() >= 3)
		{
			BOOST_CHECK_EQUAL(cbs[0][0], "AAT");
			BOOST_CHECK_EQUAL(cbs[0][1], "GAA");
			BOOST_CHECK_EQUAL(cbs[0][2], "AAA");
		}

		if (cbs[1].size() >= 3)
		{
			BOOST_CHECK_EQUAL(cbs[1][0], "TTAGGTCCA");
			BOOST_CHECK_EQUAL(cbs[1][1], "TTAGGGGCC");
			BOOST_CHECK_EQUAL(cbs[1][2], "TTAGGTCCC");
		}

		Merge::BarcodesParsing::InDropBarcodesParser parser("");
		BOOST_CHECK_THROW(parser.get_barcodes_list("/barcodes.wrong"), std::runtime_error);
	}

	BOOST_FIXTURE_TEST_CASE(testConstLengthBarcodesFile, Fixture)
	{
		Merge::BarcodesParsing::ConstLengthBarcodesParser parser(PROJ_DATA_PATH + std::string("/barcodes/hifibio"));
		parser.init();
		auto cbs = parser._barcodes;

		BOOST_REQUIRE_EQUAL(cbs.size(), 3);

		BOOST_REQUIRE_EQUAL(cbs[0].size(), 96);
		BOOST_REQUIRE_EQUAL(cbs[1].size(), cbs[0].size());
		BOOST_REQUIRE_EQUAL(cbs[2].size(), cbs[0].size());

		BOOST_CHECK_EQUAL(cbs[0][0], "GCTTGGTTGCGCTTGGTTGC");
		BOOST_CHECK_EQUAL(cbs[1][0], "ACTTGGCGCTACTTGGCGCT");
		BOOST_CHECK_EQUAL(cbs[2][0], "GCCATTGGACGCCATTGGAC");

		BOOST_CHECK_EQUAL(cbs[0][1], "GGAGTGGTTCGGAGTGGTTC");
		BOOST_CHECK_EQUAL(cbs[1][1], "AGAACGCTCCAGAACGCTCC");
		BOOST_CHECK_EQUAL(cbs[2][1], "GCCGAATCCTGCCGAATCCT");
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
		auto cell_ids_by_cb = this->container_full->cell_ids_by_cb();
		long is1 = RealBarcodesMergeStrategy::get_umigs_intersect_size(this->container_full->cell_genes(cell_ids_by_cb["AAATTAGGTCCA"]),
																		 this->container_full->cell_genes(cell_ids_by_cb["CCCTTAGGTCCA"]));

		long is2 = RealBarcodesMergeStrategy::get_umigs_intersect_size(this->container_full->cell_genes(cell_ids_by_cb["AAATTAGGTCCC"]),
																		 this->container_full->cell_genes(cell_ids_by_cb["AAATTAGGTCCG"]));

		long is3 = RealBarcodesMergeStrategy::get_umigs_intersect_size(this->container_full->cell_genes(cell_ids_by_cb["AAATTAGGTCCA"]),
																		 this->container_full->cell_genes(cell_ids_by_cb["AAATTAGGTCCC"]));

		BOOST_CHECK_EQUAL(is1, 2);
		BOOST_CHECK_EQUAL(is2, 1);
		BOOST_CHECK_EQUAL(is3, 0);
	}

	BOOST_FIXTURE_TEST_CASE(testFillDistances, Fixture)
	{
		static const std::string ar_cbs1[] = {"AAT", "AAA", "CCT"};
		std::vector<std::string> cbs1(ar_cbs1, ar_cbs1 + sizeof(ar_cbs1) / sizeof(ar_cbs1[0]));

		Merge::BarcodesParsing::InDropBarcodesParser parser("");
		Merge::BarcodesParsing::InDropBarcodesParser::barcode_parts_list_t barcodes;
		barcodes.push_back(cbs1);
		barcodes.push_back(cbs1);

		parser._barcode2_length = 3;
		parser._barcodes = barcodes;
		auto dists = parser.get_distances_to_barcode("ACTACT");

		BOOST_CHECK_EQUAL(dists[0].size(), 3);
		BOOST_CHECK_EQUAL(dists[1].size(), 3);

		BOOST_CHECK_EQUAL(dists[0][0].value, 1);
		BOOST_CHECK_EQUAL(dists[0][1].value, 1);
		BOOST_CHECK_EQUAL(dists[0][2].value, 2);
		BOOST_CHECK_EQUAL(dists[0][2].index, 1);

		BOOST_CHECK_EQUAL(dists[1][0].value, 1);
		BOOST_CHECK_EQUAL(dists[1][1].value, 1);
		BOOST_CHECK_EQUAL(dists[1][2].value, 2);
		BOOST_CHECK_EQUAL(dists[1][2].index, 1);
	}

	BOOST_FIXTURE_TEST_CASE(testRealNeighboursCbs, Fixture)
	{
		CellsDataContainer::ids_t ids = this->real_cb_strat->get_real_neighbour_cbs(*this->container_full, this->container_full->cell_ids_by_cb().at("CAATTAGGTCCG"));

		BOOST_REQUIRE_EQUAL(ids.size(), 2);

		BOOST_CHECK_EQUAL(this->container_full->cell_barcode(ids[0]), "AAATTAGGTCCA");
		BOOST_CHECK_EQUAL(this->container_full->cell_barcode(ids[1]), "AAATTAGGTCCC");

		ids = this->real_cb_strat->get_real_neighbour_cbs(*this->container_full, this->container_full->cell_ids_by_cb().at("AAATTAGGTCCC"));

		BOOST_CHECK_EQUAL(ids.size(), 1);

		if (ids.size() >= 2)
		{
			BOOST_CHECK_EQUAL(this->container_full->cell_barcode(ids[0]), "AAATTAGGTCCC");
		}
	}

	BOOST_FIXTURE_TEST_CASE(testRealNeighbours, Fixture)
	{
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 0), 0);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 1), 1);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 2), 1);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 3), 0);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 4), 0);
		BOOST_CHECK_EQUAL(this->real_cb_strat->get_merge_target(*this->container_full, 5), 0);
	}

	BOOST_FIXTURE_TEST_CASE(testMergeByRealBarcodes, Fixture)
	{
		this->container_full->merge_and_filter();

		BOOST_CHECK_EQUAL(this->container_full->cell_barcodes_raw().size(), 7);
		BOOST_CHECK_EQUAL(this->container_full->filtered_cells().size(), 2);

		BOOST_REQUIRE(this->container_full->filtered_cells().size() >= 2);

		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[0]).size(), 3);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).size(), 4);

		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[0]).at("Gene1").size(), 1);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[0]).at("Gene1").at("CAACCT"), 2);

		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene1").size(), 2);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene1").at("AAACCT"), 3);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene2").size(), 1);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene2").at("CCCCCT"), 4);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene3").size(), 2);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene3").at("ACCCCT"), 2);
		BOOST_CHECK_EQUAL(this->container_full->cell_genes(this->container_full->filtered_cells()[1]).at("Gene3").at("CCATTC"), 1);

		BOOST_CHECK(!this->container_full->is_cell_merged(0));
		BOOST_CHECK(!this->container_full->is_cell_merged(1));
		BOOST_CHECK(this->container_full->is_cell_merged(2));
		BOOST_CHECK(this->container_full->is_cell_merged(3));
		BOOST_CHECK(this->container_full->is_cell_merged(4));
		BOOST_CHECK(this->container_full->is_cell_merged(5));
		BOOST_CHECK(!this->container_full->is_cell_merged(6));

		BOOST_CHECK_EQUAL(this->container_full->excluded_cells().size(), 1);
	}

	BOOST_FIXTURE_TEST_CASE(testSplitBarcode, Fixture)
	{
		Merge::BarcodesParsing::ConstLengthBarcodesParser parser(PROJ_DATA_PATH + std::string("/barcodes/hifibio"));
		parser.init();
		auto splitted = parser.split_barcode("TAATGAGCACTAATGAGCACATACTGATGCATACTGATGCGTAACACGCTGTAACACGCT");

		BOOST_REQUIRE_EQUAL(splitted.size(), 3);
		BOOST_CHECK_EQUAL(splitted[0], "TAATGAGCACTAATGAGCAC");
		BOOST_CHECK_EQUAL(splitted[1], "ATACTGATGCATACTGATGC");
		BOOST_CHECK_EQUAL(splitted[2], "GTAACACGCTGTAACACGCT");
	}

BOOST_AUTO_TEST_SUITE_END()