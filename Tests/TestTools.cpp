#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/unordered_map.hpp>
#include <Tools/ReadParameters.h>

#include "Tools/GeneInfo.h"
#include "Tools/Logs.h"
#include "Tools/RefGenesContainer.h"
#include "Tools/UtilFunctions.h"

using namespace Tools;

struct Fixture
{
	Fixture()
	{
		init_test_logs();
	}
};

BOOST_AUTO_TEST_SUITE(TestTools)

	BOOST_FIXTURE_TEST_CASE(testGtf, Fixture)
	{
		std::string test_str = "chr1\tunknown\tCDS\t878633  878757  .       +       2       gene_id \"SAMD11\"; "
				"gene_name \"SAMD11\"; p_id \"P11277\"; transcript_id \"NM_152486\"; tss_id \"TSS28354\";";

		auto info = RefGenesContainer::parse_gtf_record(test_str);

		BOOST_CHECK_EQUAL(info.chr_name(), "chr1");
		BOOST_CHECK_EQUAL(info.id(), "SAMD11");
		BOOST_CHECK_EQUAL(info.start_pos(), 878632);
		BOOST_CHECK_EQUAL(info.end_pos(), 878757);
	}


	BOOST_FIXTURE_TEST_CASE(testEditDistance, Fixture)
	{
		BOOST_CHECK_EQUAL(Tools::edit_distance("ATTTTC", "ATTTGC"), 1);
		BOOST_CHECK_EQUAL(Tools::edit_distance("ATTTTCC", "ATTTGNC"), 1);
		BOOST_CHECK_EQUAL(Tools::edit_distance("ATTTTCC", "ATTTGNC", false), 2);
		BOOST_CHECK_EQUAL(Tools::edit_distance("ATTTTCC", "ATTTGTC"), 2);
		BOOST_CHECK_EQUAL(Tools::edit_distance("ATTTTCC", "ATTTTCC"), 0);
	}

	BOOST_FIXTURE_TEST_CASE(testReadParams, Fixture)
	{
		ReadParameters rp("@111!ATTTGC#ATATC");
		BOOST_CHECK_EQUAL(rp.read_name(), "@111");
		BOOST_CHECK_EQUAL(rp.read_name_safe(), "@111");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTGC");
		BOOST_CHECK_EQUAL(rp.umi(), "ATATC");

		rp = ReadParameters("111!ATTTG#ATAT");
		BOOST_CHECK_EQUAL(rp.read_name(), "111");
		BOOST_CHECK_EQUAL(rp.read_name_safe(), "111");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTG");
		BOOST_CHECK_EQUAL(rp.umi(), "ATAT");

		rp = ReadParameters("!ATTTGC#ATATC");
		BOOST_CHECK_EQUAL(rp.read_name(), "");
		BOOST_CHECK_EQUAL(rp.read_name_safe(), "!ATTTGC#ATATC");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTGC");
		BOOST_CHECK_EQUAL(rp.umi(), "ATATC");

		rp = ReadParameters("trash!ATTTG#ATAT");
		BOOST_CHECK_EQUAL(rp.read_name(), "trash");
		BOOST_CHECK_EQUAL(rp.read_name_safe(), "trash");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTG");
		BOOST_CHECK_EQUAL(rp.umi(), "ATAT");

		ReadParameters rp2 = ReadParameters(rp.to_monolithic_string());
		BOOST_CHECK_EQUAL(rp2.read_name(), rp.read_name());
		BOOST_CHECK_EQUAL(rp2.read_name_safe(), rp2.read_name());
		BOOST_CHECK_EQUAL(rp2.cell_barcode(), rp.cell_barcode());
		BOOST_CHECK_EQUAL(rp2.umi(), rp.umi());

		rp2 = ReadParameters(rp.to_monolithic_string("1111"));
		BOOST_CHECK_EQUAL(rp2.read_name(), "1111");
		BOOST_CHECK_EQUAL(rp2.read_name_safe(), rp2.read_name());
		BOOST_CHECK_EQUAL(rp2.cell_barcode(), rp.cell_barcode());
		BOOST_CHECK_EQUAL(rp2.umi(), rp.umi());

		BOOST_CHECK_THROW(ReadParameters("ATTTG#ATAT"), std::runtime_error);
	}

	BOOST_FIXTURE_TEST_CASE(testGeneMerge, Fixture)
	{
		IntervalsContainer<std::string> intervals;
		intervals.add_interval(0, 100, "");
		intervals.add_interval(200, 300, "");
		intervals.add_interval(400, 500, "");

		intervals.set_initialized(false);

		auto &homogenous_intervals = intervals._homogenous_intervals;
		BOOST_REQUIRE_EQUAL(homogenous_intervals.size(), 3);
		BOOST_CHECK_EQUAL(homogenous_intervals.back().start_pos(), 400);
		BOOST_CHECK_EQUAL(homogenous_intervals.back().end_pos(), 500);

		auto cur_int = Interval(90, 110);
		BOOST_CHECK_EQUAL(homogenous_intervals.front().is_intercept(cur_int), true);
		intervals.add_interval(90, 110, "", true);
		intervals.set_initialized(false);
		BOOST_CHECK_EQUAL(homogenous_intervals.front().end_pos(), 110);
		cur_int = Interval(110, 130);
		BOOST_CHECK_EQUAL(homogenous_intervals.front().is_intercept(cur_int), false);

		intervals.add_interval(150, 190, "", true);
		intervals.set_initialized(false);
		BOOST_CHECK_EQUAL(homogenous_intervals.size(), 4);

		intervals.add_interval(110, 151, "", true);
		intervals.set_initialized(false);
		BOOST_CHECK_EQUAL(homogenous_intervals.size(), 3);
		BOOST_CHECK_EQUAL(homogenous_intervals.front().end_pos(), 190);

		intervals.add_interval(190, 401, "", true);
		intervals.set_initialized(false);
		BOOST_CHECK_EQUAL(homogenous_intervals.size(), 1);
		BOOST_CHECK_EQUAL(homogenous_intervals.front().start_pos(), 0);
		BOOST_CHECK_EQUAL(homogenous_intervals.front().end_pos(), 500);
	}

	BOOST_FIXTURE_TEST_CASE(testInitGtf, Fixture)
	{
		const std::string gtf_filename = PROJ_DATA_PATH + (std::string)("/gtf/gtf_test.gtf.gz");
		RefGenesContainer genes_container(gtf_filename);
		
		BOOST_CHECK_EQUAL(genes_container._genes_intervals.size(), 3);
		auto &chr1_intervals = genes_container._genes_intervals["chr1"]._homogenous_intervals;
		BOOST_REQUIRE_EQUAL(chr1_intervals.size(), 10);

		BOOST_CHECK_EQUAL(chr1_intervals[0].start_pos(), 11873);
		BOOST_CHECK_EQUAL(chr1_intervals[0].end_pos(), 14209);
		BOOST_CHECK_EQUAL(chr1_intervals[1].start_pos(), 14361);
		BOOST_CHECK_EQUAL(chr1_intervals[1].end_pos(), 18366);
		BOOST_CHECK_EQUAL(chr1_intervals[2].start_pos(), 24320);
		BOOST_CHECK_EQUAL(chr1_intervals[2].end_pos(), 29370);
		BOOST_CHECK_EQUAL(chr1_intervals[3].start_pos(), 34610);
		BOOST_CHECK_EQUAL(chr1_intervals[3].end_pos(), 35174);
		BOOST_CHECK_EQUAL(chr1_intervals[4].start_pos(), 35276);
		BOOST_CHECK_EQUAL(chr1_intervals[4].end_pos(), 35481);
		BOOST_CHECK_EQUAL(chr1_intervals[5].start_pos(), 69090);
		BOOST_CHECK_EQUAL(chr1_intervals[5].end_pos(), 69499);

		BOOST_CHECK_EQUAL(chr1_intervals[6].start_pos(), 69499);
		BOOST_CHECK_EQUAL(chr1_intervals[6].end_pos(), 69790);
		BOOST_CHECK_EQUAL(chr1_intervals[7].start_pos(), 69790);
		BOOST_CHECK_EQUAL(chr1_intervals[7].end_pos(), 70008);
		BOOST_CHECK_EQUAL(chr1_intervals[8].start_pos(), 70008);
		BOOST_CHECK_EQUAL(chr1_intervals[8].end_pos(), 71005);
		BOOST_CHECK_EQUAL(chr1_intervals[9].start_pos(), 71005);
		BOOST_CHECK_EQUAL(chr1_intervals[9].end_pos(), 72008);

		BOOST_CHECK_EQUAL(chr1_intervals[0].base_interval_labels.size(), 1);
		BOOST_CHECK_EQUAL(chr1_intervals[3].base_interval_labels.size(), 2);
		BOOST_CHECK_EQUAL(chr1_intervals[4].base_interval_labels.size(), 2);
		BOOST_CHECK_EQUAL(chr1_intervals[6].base_interval_labels.size(), 2);
		BOOST_CHECK_EQUAL(chr1_intervals[7].base_interval_labels.size(), 3);
		BOOST_CHECK_EQUAL(chr1_intervals[8].base_interval_labels.size(), 2);
		BOOST_CHECK_EQUAL(chr1_intervals[9].base_interval_labels.size(), 1);
		BOOST_CHECK_EQUAL(chr1_intervals.size(), 10);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals["chr2"]._homogenous_intervals.size(), 5);
	}

	BOOST_FIXTURE_TEST_CASE(testParseBed, Fixture)
	{
		const std::string gtf_filename = PROJ_DATA_PATH + (std::string)("/gtf/refflat_ucsc_mm10_exons.gtf.gz");
		const std::string bed_filename = PROJ_DATA_PATH + (std::string)("/gtf/refflat_ucsc_mm10.trimmed.bed.gz");
		RefGenesContainer gtf_container(gtf_filename);
		RefGenesContainer bed_container(bed_filename);

		BOOST_REQUIRE_EQUAL(gtf_container._genes_intervals.size(), bed_container._genes_intervals.size());
		for (auto const &gtf_chr : gtf_container._genes_intervals)
		{
			auto bed_chr_it = bed_container._genes_intervals.find(gtf_chr.first);
			if (bed_chr_it == bed_container._genes_intervals.end())
			{
				std::cout << gtf_chr.first << " isn't presented in bed file" << std::endl;
				continue;
			}

			if (gtf_chr.second._homogenous_intervals.size() != bed_chr_it->second._homogenous_intervals.size())
			{
				std::cout << gtf_chr.first << " has different number of intervals in bed and gtf" << std::endl;
				continue;
			}

			int failed_num = 0;
			for (size_t i = 0; i < gtf_chr.second._homogenous_intervals.size(); ++i)
			{
				auto const &gtf_interval = gtf_chr.second._homogenous_intervals[i];
				auto const &bed_interval = bed_chr_it->second._homogenous_intervals[i];

				bool eq = (gtf_interval.base_interval_labels.size() == bed_interval.base_interval_labels.size()) &&
						  (gtf_interval.start_pos() == bed_interval.start_pos()) &&
						  (gtf_interval.end_pos() == bed_interval.end_pos());
				failed_num += !eq;
			}
			BOOST_CHECK_EQUAL(failed_num, 0);
			if (failed_num != 0)
			{
				std::cout << gtf_chr.first << ": " << 100.0 * failed_num / gtf_chr.second._homogenous_intervals.size()
				          << "% failed" << std::endl;
				continue;
			}
		}
	}

	BOOST_FIXTURE_TEST_CASE(testGeneNames, Fixture)
	{
		init_test_logs(boost::log::trivial::info);
		const std::string gtf_filename = PROJ_DATA_PATH + (std::string)("/gtf/gtf_test.gtf.gz");
		RefGenesContainer genes_container(gtf_filename);

		RefGenesContainer::gene_names_set_t res;
		res = genes_container.get_gene_info("chr1", 11874, 12627);
		BOOST_REQUIRE_EQUAL(res.size(), 1);
		BOOST_CHECK_EQUAL(*res.begin(), "DDX11L1");

		res = genes_container.get_gene_info("chr1", 17106, 17742);
		BOOST_REQUIRE_EQUAL(res.size(), 1);
		BOOST_CHECK_EQUAL(*res.begin(), "WASH7P");

		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 30000, 31000).size(), 0);

		res = genes_container.get_gene_info("chr1", 34621, 35074);
		BOOST_REQUIRE_EQUAL(res.size(), 2);
		BOOST_CHECK_EQUAL(*res.begin(), "FAM138A");
		BOOST_CHECK_EQUAL(*(++res.begin()), "FAM138F");
//		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 69791, 69793).id(), "AR4F5,BR4F5,OR4F5");
//		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 50000, 69793).id(), "AR4F5,BR4F5,OR4F5");
//		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 35277, 69793).id(), "AR4F5,BR4F5,FAM138A,FAM138F,OR4F5");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chrX", 0, 34608).size(), 0);
		BOOST_CHECK_EQUAL(*genes_container.get_gene_info("chrX", 34609, 34612).begin(), "FAM138A");
		BOOST_CHECK_EQUAL(*genes_container.get_gene_info("chrX", 34609, 35174).begin(), "FAM138A");
//		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chrX", 100000, 110000).name(), "CHRX_GENE,CHRX_GENE2");
//		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chrX", 100000, 130000).name(), "CHRX_GENE,CHRX_GENE2");

		BOOST_CHECK_THROW(genes_container.get_gene_info("chr3", 0, 100), RefGenesContainer::ChrNotFoundException);
		BOOST_CHECK_THROW(genes_container.get_gene_info("chrM", 100000, 130000), RefGenesContainer::ChrNotFoundException);
	}

	BOOST_FIXTURE_TEST_CASE(testR, Fixture)
	{
		auto R = init_r();
		R->parseEval("library(ggplot2)\n"
					 "library(grid)\n"
					 "library(gridExtra)\n"
					 "library(knitr)\n"
					 "library(parallel)\n"
					 "library(fitdistrplus)");

		R->parseEval((std::string)"source('" + PROJ_BIN_PATH + "/Functions.R')");
	}

	BOOST_FIXTURE_TEST_CASE(testInterval, Fixture)
	{
		IntervalsContainer<std::string> intervals;
		intervals.add_interval(10, 20, "i1");
		intervals.add_interval(10, 20, "i1");
		intervals.add_interval(15, 30, "i1");
		intervals.add_interval(15, 20, "i2");

		intervals.set_initialized();

		BOOST_REQUIRE_EQUAL(intervals.get_intervals(0, 11).size(), 1);
		BOOST_CHECK_EQUAL(*intervals.get_intervals(0, 11).begin(), "i1");

		BOOST_CHECK(intervals.get_intervals(0, 5).empty());

		BOOST_REQUIRE_EQUAL(intervals.get_intervals(25, 30).size(), 1);
		BOOST_CHECK_EQUAL(*intervals.get_intervals(25, 30).begin(), "i1");

		BOOST_REQUIRE_EQUAL(intervals.get_intervals(17, 20).size(), 2);
		BOOST_CHECK_EQUAL(*intervals.get_intervals(17, 20).begin(), "i1");
		BOOST_CHECK_EQUAL(*(++intervals.get_intervals(17, 20).begin()), "i2");
	}

//	BOOST_FIXTURE_TEST_CASE(testGtfPerformance, Fixture) //Uncomment to print performance
//	{
//		init_test_logs(boost::log::trivial::info);
//
//		time_t t0 = clock();
//		std::cout << "Start init " << std::endl;
//		RefGenesContainer genes_container(PROJ_DATA_PATH + (std::string) "/gtf/gencode.v19.annotation.gtf.gz");
//		double t_len = (clock() - t0) / (CLOCKS_PER_SEC / 1000.0);
//		std::cout << "End init: " << t_len << "ms" << std::endl;
//
//		t_len = 0;
//		srand(time(0));
//		std::cout << "Start 100000 search" << std::endl;
//		for (int i = 0; i < 100000; ++i)
//		{
//			int sp = rand() % 1000000;
//			t0 = clock();
//			genes_container.get_gene_info("chr10", sp, sp + rand() % 1000000);
//			t_len += (clock() - t0) / (CLOCKS_PER_SEC / 1000.0);
//		}
//		std::cout << "End 100000 search: " << t_len << "ms" << std::endl;
//	}

BOOST_AUTO_TEST_SUITE_END()
