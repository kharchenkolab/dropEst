#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <boost/unordered_map.hpp>

#include "Tools/RefGenesContainer.h"
#include "Tools/Logs.h"

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
		BOOST_CHECK_EQUAL(info.start_pos(), 878633);
		BOOST_CHECK_EQUAL(info.end_pos(), 878757);
	}

	BOOST_FIXTURE_TEST_CASE(testGeneMerge, Fixture)
	{
		std::list<RefGenesContainer::GeneInfo> genes;
		RefGenesContainer::GeneInfo inf;
		inf._start_pos = 0; inf._end_pos = 100;
		genes.push_back(inf);
		inf._start_pos = 200; inf._end_pos = 300;
		genes.push_back(inf);
		inf._start_pos = 400; inf._end_pos = 500;
		genes.push_back(inf);

		RefGenesContainer::add_gene(inf, genes);
		BOOST_CHECK_EQUAL(genes.size(), 3);
		BOOST_CHECK_EQUAL(genes.back()._start_pos, 400);
		BOOST_CHECK_EQUAL(genes.back()._end_pos, 500);

		inf._start_pos = 90; inf._end_pos = 110;
		BOOST_CHECK_EQUAL(genes.front().is_intercept(inf), true);
		genes.front().merge(inf);
		BOOST_CHECK_EQUAL(genes.front()._end_pos, 110);

		inf._start_pos = 150; inf._end_pos = 190;
		RefGenesContainer::add_gene(inf, genes);
		BOOST_CHECK_EQUAL(genes.size(), 4);

		inf._start_pos = 110; inf._end_pos = 150;
		RefGenesContainer::add_gene(inf, genes);
		BOOST_CHECK_EQUAL(genes.size(), 3);
		BOOST_CHECK_EQUAL(genes.front()._end_pos, 190);

		inf._start_pos = 190; inf._end_pos = 400;
		RefGenesContainer::add_gene(inf, genes);
		BOOST_CHECK_EQUAL(genes.size(), 1);
		BOOST_CHECK_EQUAL(genes.front()._start_pos, 0);
		BOOST_CHECK_EQUAL(genes.front()._end_pos, 500);
	}

	BOOST_FIXTURE_TEST_CASE(testInitGtf, Fixture)
	{
		RefGenesContainer genes_container;
		const std::string gtf_filename = PROJ_DATA_PATH + (std::string)("/gtf_test.gtf");
		genes_container.init_from_gtf(gtf_filename);

		BOOST_CHECK_EQUAL(genes_container._genes_intervals.size(), 3);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][0].start_pos, 11874);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][0].end_pos, 14209);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][1].start_pos, 14362);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][1].end_pos, 18366);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][2].start_pos, 24321);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][2].end_pos, 29370);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][3].start_pos, 34611);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][3].end_pos, 35174);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][4].start_pos, 35277);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][4].end_pos, 35481);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][5].start_pos, 69091);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][5].end_pos, 69500);

		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][6].start_pos, 69500);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][6].end_pos, 69791);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][7].start_pos, 69791);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][7].end_pos, 70008);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][8].start_pos, 70008);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][8].end_pos, 71005);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][9].start_pos, 71005);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][9].end_pos, 72008);

		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][0].genes.size(), 1);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][3].genes.size(), 2);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][4].genes.size(), 2);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][6].genes.size(), 2);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][7].genes.size(), 3);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][8].genes.size(), 2);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1][9].genes.size(), 1);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[1].size(), 10);
		BOOST_CHECK_EQUAL(genes_container._genes_intervals[2].size(), 2);
	}

	BOOST_FIXTURE_TEST_CASE(testGenesNames, Fixture)
	{
		RefGenesContainer genes_container;
		const std::string gtf_filename = PROJ_DATA_PATH + (std::string)("/gtf_test.gtf");
		genes_container.init_from_gtf(gtf_filename);

		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 11874, 12627).id(), "DDX11L1");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 17106, 17742).id(), "WASH7P");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 17106, 24738).id(), "WASH7P");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 30000, 31000).id(), "");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 34621, 35074).id(), "FAM138A,FAM138F");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 69791, 69793).id(), "AR4F5,BR4F5,OR4F5");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 50000, 69793).id(), "AR4F5,BR4F5,OR4F5");
		BOOST_CHECK_EQUAL(genes_container.get_gene_info("chr1", 35277, 69793).id(), "AR4F5,BR4F5,FAM138A,FAM138F,OR4F5");
		BOOST_CHECK_NO_THROW(genes_container.get_gene_info("chr1", 35277, 105000));
	}

//	BOOST_FIXTURE_TEST_CASE(testGtfPerformance, Fixture) //Uncomment to print performance
//	{
//		init_test_logs(boost::log::trivial::info);
//		RefGenesContainer genes_container;
//
//		time_t t0 = clock();
//		std::cout << "Start init" << std::endl;
//		genes_container.init_from_gtf("/home/victor/Study/InDrop/Data/genes.gtf");
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
//			genes_container.get_gene_info("chr3", sp, sp + rand() % 1000000);
//			t_len += (clock() - t0) / (CLOCKS_PER_SEC / 1000.0);
//		}
//		std::cout << "End 100000 search: " << t_len << "ms" << std::endl;
//	}

BOOST_AUTO_TEST_SUITE_END()