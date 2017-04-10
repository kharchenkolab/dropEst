#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Estimation/CellsDataContainer.h"

#include <RInside.h>
#include <Estimation/BamProcessing/ReadsParamsParser.h>
#include <Estimation/BamProcessing/BamController.h>
#include <Estimation/BamProcessing/BamProcessor.h>
#include <Estimation/Merge/BarcodesParsing/InDropBarcodesParser.h>
#include <Estimation/Merge/BarcodesParsing/ConstLengthBarcodesParser.h>
#include <Estimation/Merge/MergeStrategyFactory.h>
#include <Estimation/Merge/SimpleMergeStrategy.h>
#include <Estimation/Merge/RealBarcodesMergeStrategy.h>
#include <Estimation/MergeUMIs/MergeUMIsStrategySimple.h>

#include <Tools/Logs.h>
#include <Tools/UtilFunctions.h>

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
		this->umi_merge_strat = std::make_shared<MergeUMIs::MergeUMIsStrategySimple>(1);
		this->container_full = std::make_shared<CellsDataContainer>(this->real_cb_strat, this->umi_merge_strat, 1);

		Tools::init_test_logs(boost::log::trivial::info);
		this->container_full->add_record("AAATTAGGTCCA", "AAACCT", "Gene1"); //0, real
		this->container_full->add_record("AAATTAGGTCCA", "CCCCCT", "Gene2");
		this->container_full->add_record("AAATTAGGTCCA", "ACCCCT", "Gene3");
		this->container_full->add_record("AAATTAGGTCCA", "ACCCCT", "Gene4");

		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene1"); //1, real
		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene10");
		this->container_full->add_record("AAATTAGGTCCC", "CAACCT", "Gene20");

		this->container_full->add_record("AAATTAGGTCCG", "CAACCT", "Gene1"); //2, false

		this->container_full->add_record("AAATTAGGTCGG", "AAACCT", "Gene1"); //3, false
		this->container_full->add_record("AAATTAGGTCGG", "CCCCCT", "Gene2");

		this->container_full->add_record("CCCTTAGGTCCA", "CCATTC", "Gene3"); //4, false
		this->container_full->add_record("CCCTTAGGTCCA", "CCCCCT", "Gene2");
		this->container_full->add_record("CCCTTAGGTCCA", "ACCCCT", "Gene3");

		this->container_full->add_record("CAATTAGGTCCG", "CAACCT", "Gene1"); //5, false
		this->container_full->add_record("CAATTAGGTCCG", "AAACCT", "Gene1");
		this->container_full->add_record("CAATTAGGTCCG", "CCCCCT", "Gene2");

		this->container_full->add_record("AAAAAAAAAAAA", "CCCCCT", "Gene2"); //6, false, excluded
		this->container_full->set_initialized();
	}

	std::shared_ptr<Merge::RealBarcodesMergeStrategy> real_cb_strat;
	std::shared_ptr<MergeUMIs::MergeUMIsStrategySimple> umi_merge_strat;
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
		std::unordered_map<std::string, int> map;
		auto &t = ++map["chr1"];
		BOOST_CHECK_EQUAL(map["chr1"], 1);
		map.emplace(std::make_pair("chr1", 0));
		BOOST_CHECK_EQUAL(map["chr1"], 1);
		t = 10;
		BOOST_CHECK_EQUAL(map["chr1"], 10);

		auto &t2 = map["chr2"];
		BOOST_CHECK_EQUAL(map["chr2"], 0);
		++t2;
		BOOST_CHECK_EQUAL(map["chr2"], 1);
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

	BOOST_FIXTURE_TEST_CASE(testStat, Fixture)
	{
		Stats stats;
		stats.set(Stats::StrStrStatType::MERGE_EDIT_DISTANCE_BY_CELL, "AAA", "BBB", 10);
		stats.set(Stats::StrStrStatType::MERGE_EDIT_DISTANCE_BY_CELL, "CCC", "BBB", 20);

		auto map = stats.get_raw(Stats::StrStrStatType::MERGE_EDIT_DISTANCE_BY_CELL);
		BOOST_REQUIRE_EQUAL(map.size(), 2);

		BOOST_CHECK_EQUAL(map["AAA"]["BBB"], 10);
		BOOST_CHECK_EQUAL(map["CCC"]["BBB"], 20);
	}

	BOOST_FIXTURE_TEST_CASE(testGeneMatchLevel, Fixture)
	{
		using namespace BamProcessing;
		BamTools::BamAlignment align;
		align.Position = 34610;
		align.Length = 10;
		align.CigarData.push_back(BamTools::CigarOp('M', 10));
		ReadsParamsParser parser0(PROJ_DATA_PATH + (std::string)"/gtf/gtf_test.gtf.gz", ReadsParamsParser::ANY);
		ReadsParamsParser parser1(PROJ_DATA_PATH + (std::string)"/gtf/gtf_test.gtf.gz", ReadsParamsParser::ONE);
		ReadsParamsParser parser2(PROJ_DATA_PATH + (std::string)"/gtf/gtf_test.gtf.gz", ReadsParamsParser::BOTH);

		BOOST_CHECK_EQUAL(parser0.get_gene("chrX", align), "FAM138A");
		BOOST_CHECK_EQUAL(parser1.get_gene("chrX", align), "");
		BOOST_CHECK_EQUAL(parser2.get_gene("chrX", align), "FAM138A");

		align.Position = 34600;
		BOOST_CHECK_EQUAL(parser0.get_gene("chrX", align), "FAM138A");
		BOOST_CHECK_EQUAL(parser1.get_gene("chrX", align), "FAM138A");
		BOOST_CHECK_THROW(parser2.get_gene("chrX", align), ReadsParamsParser::MoleculeHasIntons);
	}

	BOOST_FIXTURE_TEST_CASE(testUmiExclusion, Fixture)
	{
		CellsDataContainer container(this->real_cb_strat, this->umi_merge_strat, 10);
		container.add_record("AAATTAGGTCCA", "AAACCT", "Gene1");
		container.add_record("AAATTAGGTCCA", "CCCCCT", "Gene2");
		container.add_record("AAATTAGGTCCA", "ACCCCT", "Gene3");
		container.add_record("AAATTAGGTCCA", "ACCCCT", "Gene4");

		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene4").at("ACCCCT"), 1);

		container.add_record("AAATTAGGTCCA", "TTTTTT", "Gene3", true);
		container.add_record("AAATTAGGTCCA", "ACCCCT", "Gene4", true);
		BOOST_CHECK(container.cell_genes(0).at("Gene3").at("TTTTTT") < 0);
		BOOST_CHECK(container.cell_genes(0).at("Gene4").at("ACCCCT") < 0);

		container.set_initialized();
		container.merge_and_filter();

		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene3").at("ACCCCT"), 1);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene3").size(), 1);

		BOOST_CHECK_THROW(container.cell_genes(0).at("Gene4").at("ACCCCT"), std::out_of_range);
		BOOST_CHECK_THROW(container.cell_genes(0).at("Gene4"), std::out_of_range);
	}

	BOOST_FIXTURE_TEST_CASE(testGeneMatchLevelUmiExclusion, Fixture)
	{
		using namespace BamProcessing;
		CellsDataContainer container(this->real_cb_strat, this->umi_merge_strat, 10);
		auto parser = std::make_shared<ReadsParamsParser>(PROJ_DATA_PATH + (std::string)"/gtf/gtf_test.gtf.gz", ReadsParamsParser::BOTH);
		std::shared_ptr<BamProcessorAbstract> processor(new BamProcessor(container, false));
		std::unordered_set<std::string> unexpected_chromosomes;

		BamTools::BamAlignment align;
		align.Name = "152228477!TGAGTTCTGTTACTGCATC#ATGGGC";
		align.Position = 34610;
		align.Length = 10;
		align.CigarData.push_back(BamTools::CigarOp('M', 10));

		BamController::process_alignment(parser, processor, unexpected_chromosomes, "chrX", align);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("FAM138A").at("ATGGGC"), 1);

		align.Position = 34600;
		BamController::process_alignment(parser, processor, unexpected_chromosomes, "chrX", align);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("FAM138A").at("ATGGGC"), CellsDataContainer::UMI_EXCLUDED);

		align.Position = 34610;
		BamController::process_alignment(parser, processor, unexpected_chromosomes, "chrX", align);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("FAM138A").at("ATGGGC"), CellsDataContainer::UMI_EXCLUDED);

		align.Name = "152228477!TGAGTTCTGTTACTGCATC#ATTTTC";
		align.Position = 34600;
		BamController::process_alignment(parser, processor, unexpected_chromosomes, "chrX", align);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("FAM138A").at("ATTTTC"), CellsDataContainer::UMI_EXCLUDED);

		container.set_initialized();
		container.merge_and_filter();

		BOOST_CHECK_THROW(container.cell_genes(0).at("FAM138A"), std::out_of_range);
	}

	BOOST_FIXTURE_TEST_CASE(testUMIMerge, Fixture)
	{
		CellsDataContainer container(this->real_cb_strat, this->umi_merge_strat, 0);

		container.add_record("AAATTAGGTCCA", "AAACCT", "Gene1");
		container.add_record("AAATTAGGTCCA", "CCCCCT", "Gene1");
		container.add_record("AAATTAGGTCCA", "AAATTN", "Gene1");
		container.add_record("AAATTAGGTCCA", "ACCCCT", "Gene1");

		CellsDataContainer::s_s_hash_t merge_targets;
		merge_targets["AAACCT"] = "CCCCCT";
		merge_targets["AAATTN"] = "GGGGGG";
		merge_targets["ACCCCT"] = "ACCCCT";

		container.merge_umis(0, "Gene1", merge_targets);

		BOOST_REQUIRE_EQUAL(container.cell_genes(0).at("Gene1").size(), 3);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene1").at("CCCCCT"), 2);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene1").at("GGGGGG"), 1);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene1").at("ACCCCT"), 1);
	}

	BOOST_FIXTURE_TEST_CASE(testFillWrongUmis, Fixture)
	{
		using namespace MergeUMIs;
		MergeUMIsStrategySimple::s_vec_t wrong_umis;
		wrong_umis.push_back("AAANTTT");
		wrong_umis.push_back("AAANCTT");
		wrong_umis.push_back("NNNNNNN");
		auto merge_targets = this->umi_merge_strat->fill_wrong_umis(wrong_umis);
		for (auto const &umi : merge_targets)
		{
			BOOST_CHECK_NE(umi.first, umi.second);
			BOOST_CHECK_EQUAL(Tools::hamming_distance(umi.first, umi.second), 0);
			BOOST_CHECK_EQUAL(umi.second.find('N'), std::string::npos);
		}
	}

	BOOST_FIXTURE_TEST_CASE(testRemoveSimilarWrongUmis, Fixture)
	{
		using namespace MergeUMIs;
		MergeUMIsStrategySimple::s_vec_t wrong_umis;
		wrong_umis.push_back("AAANTTT");
		wrong_umis.push_back("AACTNTT");
		wrong_umis.push_back("AACTCNT");
//		wrong_umis.push_back("NNNNNNN"); // These tests will fail, but we don't really care
//
//		wrong_umis.push_back("AANNNCT");
//		wrong_umis.push_back("AGNNNCT");
//		wrong_umis.push_back("ANATTCT");
		this->umi_merge_strat->remove_similar_wrong_umis(wrong_umis);
		BOOST_CHECK_EQUAL(wrong_umis.size(), 2);
		BOOST_CHECK_EQUAL(wrong_umis[0], "AAANTTT");
		BOOST_CHECK_EQUAL(wrong_umis[1], "AACTCNT");
	}

	BOOST_FIXTURE_TEST_CASE(testUMIMergeStrategy, Fixture)
	{
		using namespace MergeUMIs;
		CellsDataContainer container(this->real_cb_strat, this->umi_merge_strat, 0);

		container.add_record("AAATTAGGTCCA", "AAACCT", "Gene1");
		container.add_record("AAATTAGGTCCA", "CCCCCT", "Gene1");
		container.add_record("AAATTAGGTCCA", "AAACCN", "Gene1");
		container.add_record("AAATTAGGTCCA", "ACCCCT", "Gene1");

		container.add_record("AAATTAGGTCCA", "TTTTTT", "Gene2");
		container.add_record("AAATTAGGTCCA", "TTTNNG", "Gene2");
		container.add_record("AAATTAGGTCCA", "TTGNNG", "Gene2");
		container.add_record("AAATTAGGTCCA", "ACCCCT", "Gene2");
		container.add_record("AAATTAGGTCCA", "NNNNNN", "Gene2");

		container.set_initialized();

		this->umi_merge_strat->merge(container);
		BOOST_REQUIRE_EQUAL(container.cell_genes(0).at("Gene1").size(), 3);
		BOOST_REQUIRE_EQUAL(container.cell_genes(0).at("Gene2").size(), 3);

		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene1").at("AAACCT"), 2);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene1").at("CCCCCT"), 1);
		BOOST_CHECK_EQUAL(container.cell_genes(0).at("Gene1").at("ACCCCT"), 1);

		BOOST_CHECK_NO_THROW(container.cell_genes(0).at("Gene2").at("TTTTTT"));
		BOOST_CHECK_NO_THROW(container.cell_genes(0).at("Gene2").at("ACCCCT"));

		for (auto const &umi :container.cell_genes(0).at("Gene2"))
		{
			BOOST_CHECK_EQUAL(umi.first.find('N'), std::string::npos);
		}
	}
BOOST_AUTO_TEST_SUITE_END()