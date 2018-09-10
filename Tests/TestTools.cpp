#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <iostream>
#include <fstream>
#include <boost/test/included/unit_test.hpp>
#include <boost/unordered_map.hpp>
#include <Tools/ReadParameters.h>
#include <Tools/GeneAnnotation/RefGenesContainer.h>

#include "Tools/GeneAnnotation/GtfRecord.h"
#include "Tools/Logs.h"
#include "Tools/GeneAnnotation/RefGenesContainer.h"
#include "Tools/UtilFunctions.h"

using namespace Tools;
using namespace Tools::GeneAnnotation;

struct Fixture
{
	Fixture()
		: test_gtf_name(PROJ_DATA_PATH + (std::string)("/gtf/gtf_test.gtf.gz"))
	{
		init_test_logs();
	}

	std::string test_gtf_name;
};

BOOST_AUTO_TEST_SUITE(TestTools)

	BOOST_FIXTURE_TEST_CASE(testGtf, Fixture)
	{
		std::string test_str = "chr1\tunknown\texon\t878633  878757  .       +       2       gene_id \"SAMD11\"; "
				"gene_name \"SAMD11\"; p_id \"P11277\"; transcript_id \"NM_152486\"; tss_id \"TSS28354\";";

		RefGenesContainer container(this->test_gtf_name);
		auto info = container.parse_gtf_record(test_str);

		BOOST_CHECK_EQUAL(info.chr_name(), "chr1");
		BOOST_CHECK_EQUAL(info.gene_id(), "SAMD11");
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
		ReadParameters rp(ReadParameters::parse_encoded_id("@111!ATTTGC#ATATC"));
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTGC");
		BOOST_CHECK_EQUAL(rp.umi(), "ATATC");

		rp = ReadParameters::parse_encoded_id("111!ATTTG#ATAT");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTG");
		BOOST_CHECK_EQUAL(rp.umi(), "ATAT");

		rp = ReadParameters::parse_encoded_id("!ATTTGC#ATATC");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTGC");
		BOOST_CHECK_EQUAL(rp.umi(), "ATATC");

		rp = ReadParameters::parse_encoded_id("trash!ATTTG#ATAT");
		BOOST_CHECK_EQUAL(rp.cell_barcode(), "ATTTG");
		BOOST_CHECK_EQUAL(rp.umi(), "ATAT");

		rp = ReadParameters("AAATTTTATA", "TTGG", "QUALCB", "CCC");
		auto rp_info = ReadParameters::parse_from_string(rp.to_string("ID"));
		BOOST_CHECK_EQUAL(rp_info.first, "ID");
		BOOST_CHECK_EQUAL(rp_info.second.cell_barcode(), rp.cell_barcode());
		BOOST_CHECK_EQUAL(rp_info.second.umi(), rp.umi());
		BOOST_CHECK_EQUAL(rp_info.second.cell_barcode_quality(), rp.cell_barcode_quality());
		BOOST_CHECK_EQUAL(rp_info.second.umi_quality(), rp.umi_quality());

		auto rp2 = ReadParameters::parse_encoded_id(rp.encoded_id("1111"));
		BOOST_CHECK_EQUAL(rp2.cell_barcode(), rp.cell_barcode());
		BOOST_CHECK_EQUAL(rp2.umi(), rp.umi());

		BOOST_CHECK_THROW(ReadParameters::parse_encoded_id("ATTTG#ATAT"), std::runtime_error);
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
		RefGenesContainer genes_container(this->test_gtf_name);

		BOOST_CHECK_EQUAL(genes_container._transcript_intervals.size(), 3);
		auto &chr1_transcripts = genes_container._transcript_intervals.at("chr1")._homogenous_intervals;
		auto &chr1_exon_intervals = genes_container._exons_by_transcripts.at("chr1");

		BOOST_REQUIRE_EQUAL(chr1_transcripts.size(), 8);

		BOOST_CHECK_EQUAL(chr1_transcripts[0].start_pos(), 11873);
		BOOST_CHECK_EQUAL(chr1_transcripts[0].end_pos(), 14209);
		BOOST_CHECK_EQUAL(chr1_transcripts[0].base_interval_labels.size(), 1);
		BOOST_CHECK_EQUAL(*chr1_transcripts[0].base_interval_labels.begin(), "NR_046018");
		BOOST_REQUIRE_EQUAL(chr1_exon_intervals.at("NR_046018")._homogenous_intervals.size(), 1);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_046018")._homogenous_intervals[0].start_pos(), 11873);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_046018")._homogenous_intervals[0].end_pos(), 14209);

		BOOST_CHECK_EQUAL(chr1_transcripts[1].start_pos(), 14361);
		BOOST_CHECK_EQUAL(chr1_transcripts[1].end_pos(), 29370);
		BOOST_CHECK_EQUAL(chr1_transcripts[1].base_interval_labels.size(), 1);
		BOOST_CHECK_EQUAL(*chr1_transcripts[1].base_interval_labels.begin(), "NR_024540");
		BOOST_REQUIRE_EQUAL(chr1_exon_intervals.at("NR_024540")._homogenous_intervals.size(), 2);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_024540")._homogenous_intervals[0].start_pos(), 14361);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_024540")._homogenous_intervals[0].end_pos(), 18366);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_024540")._homogenous_intervals[1].start_pos(), 24320);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_024540")._homogenous_intervals[1].end_pos(), 29370);

		BOOST_CHECK_EQUAL(chr1_transcripts[2].start_pos(), 34610);
		BOOST_CHECK_EQUAL(chr1_transcripts[2].end_pos(), 35481);
		BOOST_CHECK_EQUAL(chr1_transcripts[2].base_interval_labels.size(), 2);
		BOOST_CHECK_EQUAL(*chr1_transcripts[2].base_interval_labels.begin(), "NR_026818_1");
		BOOST_CHECK_EQUAL(*(++chr1_transcripts[2].base_interval_labels.begin()), "NR_026820_1");
		BOOST_REQUIRE_EQUAL(chr1_exon_intervals.at("NR_026818_1")._homogenous_intervals.size(), 2);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_026818_1")._homogenous_intervals[0].start_pos(), 34610);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_026818_1")._homogenous_intervals[0].end_pos(), 35174);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_026818_1")._homogenous_intervals[1].start_pos(), 35276);
		BOOST_CHECK_EQUAL(chr1_exon_intervals.at("NR_026818_1")._homogenous_intervals[1].end_pos(), 35481);

		BOOST_CHECK_EQUAL(chr1_transcripts[3].start_pos(), 69090);
		BOOST_CHECK_EQUAL(chr1_transcripts[3].end_pos(), 69499);
		BOOST_CHECK_EQUAL(chr1_transcripts[3].base_interval_labels.size(), 1); //ORF45
		BOOST_CHECK_EQUAL(chr1_transcripts[4].start_pos(), 69499);
		BOOST_CHECK_EQUAL(chr1_transcripts[4].end_pos(), 69790);
		BOOST_CHECK_EQUAL(chr1_transcripts[4].base_interval_labels.size(), 2);  //ORF45, ARF45
		BOOST_CHECK_EQUAL(chr1_transcripts[5].start_pos(), 69790);
		BOOST_CHECK_EQUAL(chr1_transcripts[5].end_pos(), 70008);
		BOOST_CHECK_EQUAL(chr1_transcripts[5].base_interval_labels.size(), 3);  //ORF45, ARF45, BRF45
		BOOST_CHECK_EQUAL(chr1_transcripts[6].start_pos(), 70008);
		BOOST_CHECK_EQUAL(chr1_transcripts[6].end_pos(), 71005);
		BOOST_CHECK_EQUAL(chr1_transcripts[6].base_interval_labels.size(), 2);  //ARF45, BRF45
		BOOST_CHECK_EQUAL(chr1_transcripts[7].start_pos(), 71005);
		BOOST_CHECK_EQUAL(chr1_transcripts[7].end_pos(), 72008);
		BOOST_CHECK_EQUAL(chr1_transcripts[7].base_interval_labels.size(), 1);  //BRF45

		BOOST_CHECK_EQUAL(genes_container._transcript_intervals.at("chr2")._homogenous_intervals.size(), 5);
	}

	BOOST_FIXTURE_TEST_CASE(testParseBed, Fixture)
	{
		const std::string gtf_filename = PROJ_DATA_PATH + (std::string)("/gtf/refflat_ucsc_mm10_exons.gtf.gz");
		const std::string bed_filename = PROJ_DATA_PATH + (std::string)("/gtf/refflat_ucsc_mm10.trimmed.bed.gz");
		RefGenesContainer gtf_container(gtf_filename);
		RefGenesContainer bed_container(bed_filename);

		srand(10);
		int total_gene_num = 0;
		for (int i = 0; i < 1000000; ++i)
		{
			int start_pos = rand() % 7000000 + 3000000;

			auto r_gtf = gtf_container.get_gene_info("chr1", start_pos, start_pos + 1);
			auto r_bed = bed_container.get_gene_info("chr1", start_pos, start_pos + 1);

			RefGenesContainer::query_results_t r_bed_exons;

			int gene_num = 0;
			for (auto const &gene_rec : r_bed)
			{
				if (gene_rec.type == GtfRecord::EXON)
				{
					r_bed_exons.insert(gene_rec);
				}
			}

			for (auto const &gene_rec : r_gtf)
			{
				if (gene_rec.type != GtfRecord::EXON)
					continue;

				gene_num++;
				total_gene_num++;
				BOOST_CHECK(r_bed_exons.find(RefGenesContainer::QueryResult(gene_rec.gene_name)) != r_bed.end());
			}

			BOOST_CHECK_EQUAL(gene_num, r_bed_exons.size());
		}
		std::cout << "Total genes: " << total_gene_num << std::endl;
	}

	BOOST_FIXTURE_TEST_CASE(testR, Fixture)
	{
		auto R = init_r();
		R->parseEval("library(ggplot2)\n"
					 "library(grid)\n"
					 "library(gridExtra)\n"
					 "library(knitr)\n"
					 "library(parallel)");

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

	BOOST_FIXTURE_TEST_CASE(testGenesWithIntrons, Fixture)
	{
		RefGenesContainer container(this->test_gtf_name);
		auto record = container.get_gene_info("chr1", 20000, 20010);
		BOOST_REQUIRE_EQUAL(record.size(), 1);
		BOOST_CHECK_EQUAL(record.begin()->gene_name, "WASH7P");
		BOOST_CHECK_EQUAL(record.begin()->type, GtfRecord::INTRON);

		record = container.get_gene_info("chr1", 24750, 24760);
		BOOST_REQUIRE_EQUAL(record.size(), 1);
		BOOST_CHECK_EQUAL(record.begin()->gene_name, "WASH7P");
		BOOST_CHECK_EQUAL(record.begin()->type, GtfRecord::EXON);

		record = container.get_gene_info("chr1", 10, 20);
		BOOST_CHECK_EQUAL(record.size(), 0);

//		record = container.get_gene_info("chr1", 23000, 24750); // TODO: uncomment it after implementation of CIGAR parsing
//		BOOST_REQUIRE_EQUAL(record.size(), 2);
//		BOOST_CHECK_EQUAL(record[0].gene_name, "WASH7P");
//		BOOST_CHECK_EQUAL(record[0].type, GtfRecord::EXON);
//		BOOST_CHECK_EQUAL(record[1].gene_name, "WASH7P");
//		BOOST_CHECK_EQUAL(record[1].type, GtfRecord::INTRON);
	}

	BOOST_FIXTURE_TEST_CASE(testRelativePaths, Fixture)
	{
		std::string fname = "../data/barcodes/indrop_v1_2";
		std::string source_fname = "/home/user/InDrop/dropEst/configs/indrop_v1_2.xml";
		std::string abs_fname = "/home/user/InDrop/dropEst/configs/../data/barcodes/indrop_v1_2";

		BOOST_CHECK_EQUAL(expand_relative_path("", ""), "");
		BOOST_CHECK_EQUAL(expand_relative_path(source_fname, ""), "");
		BOOST_CHECK_EQUAL(expand_relative_path(source_fname, fname), abs_fname);
		BOOST_CHECK_EQUAL(expand_relative_path(source_fname, abs_fname), abs_fname);
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
