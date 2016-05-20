#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <boost/test/unit_test.hpp>

#include "TagsSearch/SpacerFinder.h"
#include "TagsSearch/TagsFinder.h"
#include "Tools/Logs.h"
#include "Tools/ReadParameters.h"

#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace TagsSearch;

struct Fixture
{
	Fixture()
	{
		Tools::init_test_logs();
		std::stringstream config;
		config << "<config>\n"
				"    <SpacerSearch>\n"
				"        <spacer>GAGTGATTGCTTGTGACGCCTT</spacer>\n"
				"        <spacer2>GAGTGATTGCTTGTGCCGCCTT</spacer2>\n"
				"        <max_spacer_edit_distance>3</max_spacer_edit_distance>\n"
				"        <spacer_prefix_length>5</spacer_prefix_length>\n"
				"        <spacer_min_pos>8</spacer_min_pos>\n"
				"        <spacer_max_pos>11</spacer_max_pos>\n"
				"        <barcode_length>8</barcode_length>\n"
				"        <umi_length>6</umi_length>\n"
				"        <r1_rc_length>8</r1_rc_length>\n"
				"    </SpacerSearch>\n"
				"    <TailTrimming>\n"
				"        <min_align_length>10</min_align_length>\n"
				"        <max_reads>10000000</max_reads>\n"
				"        <poly_a_tail>AAAAAAAA</poly_a_tail>\n"
				"    </TailTrimming>\n"
				"</config>";

		boost::property_tree::ptree pt;
		read_xml(config, pt);

		this->spacer_finder = SpacerFinder(pt.get_child("config.SpacerSearch"));
		this->tags_finder = TagsFinder(this->spacer_finder, pt.get_child("config.TailTrimming"));
	}

	SpacerFinder spacer_finder;
	TagsFinder tags_finder;
};

BOOST_AUTO_TEST_SUITE(TestTagsSearch)

	BOOST_FIXTURE_TEST_CASE(test1, Fixture)
	{
		std::string r1_line2 = "TTCGGTTCGGAGTGATTGCTTGTGACGCCTTCTTCGATTCGCCATTTTTTTTTTT";
		std::string r2_line2 = "TTGTTTCGCCCGGTTTTCTGTTTTCAGTAAAGTCTCGTTACGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

		SpacerFinder::len_t spacer_pos = spacer_finder.find_spacer(r1_line2);
		TagsFinder::len_t r2_trim = tags_finder.get_trim_position(spacer_pos, r1_line2, r2_line2);

		BOOST_CHECK_EQUAL(spacer_pos, 9);
		BOOST_CHECK_EQUAL(r2_trim, 44);
	}

	BOOST_FIXTURE_TEST_CASE(test2, Fixture)
	{
		std::string r1_seq = "TAGTTTCGGAGTGTTTGCTTGTGACGCCTTACCTTGCCCGCGACTTTTTTTTTTT";
		std::string r2_seq = "TCTTCCACTAATAGTTATGTCATCCCTCTTATTAATCATCATCCTAGCCCTAAGTCTGGCCTATGAGTCACTACAAAAAGGATTAGACTGAACCG";
		std::string r2_description = "+";
		std::string r2_quality_str = "1>111@1@111@33AA3BAA33DE1AA0FF3DA33AB3AF3D2A12110AB000DFGD01F10A121A11A2BFB110/AA0ABG111A111BF>";
		Tools::ReadParameters res = tags_finder.fill_parameters(0, r1_seq, "r2_id", r2_seq, r2_description, r2_quality_str);
		BOOST_CHECK_EQUAL(res.is_empty(), false);
	}

	BOOST_FIXTURE_TEST_CASE(test3, Fixture)
	{
		std::string r1_seq = "TGACCATTACTGAGTGATTGCTTGTGACGCCTTAAGCGTACAGATTATTTT";
		std::string r2_seq = "GACTGGTTGAAATTGATGATTGACATTAATAATGA";
		std::string r2_description = "+";
		std::string r2_quality_str = "1>111@1@111@33AA3BAA33DE1AA0FF3DA33AB3AF3D2A12110AB000DFGD01F10A121A11A2BFB110/AA0ABG111A111BF>";
		Tools::ReadParameters res = tags_finder.fill_parameters(0, r1_seq, "r2_id", r2_seq, r2_description, r2_quality_str);
		BOOST_CHECK_EQUAL(res.is_empty(), false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestSpacerFinder)

	BOOST_FIXTURE_TEST_CASE(test1, Fixture)
	{
		std::string r1_line2 = "TAGTCTAGGAGTGATTGCTTGTGACGCCTTTCATCCTTATAATATTTTTTTTTTT";
		SpacerFinder::len_t spacer_pos = spacer_finder.find_spacer(r1_line2);
		BOOST_CHECK_EQUAL(spacer_pos, 8);
		BOOST_CHECK_EQUAL(spacer_finder.parse_cell_barcode(r1_line2, spacer_pos), "TAGTCTAGTCATCCTT");
		BOOST_CHECK_EQUAL(spacer_finder.parse_umi_barcode(r1_line2, spacer_pos), "ATAATA");
	}

	BOOST_FIXTURE_TEST_CASE(test2, Fixture)
	{
		std::string r1_line2 = "TAGTTTCGGAGTGTTTGCTTGTGACGCCTTACCTTGCCCGCGACTTTTTTTTTTT";
		SpacerFinder::len_t spacer_pos = spacer_finder.find_spacer(r1_line2);
		BOOST_CHECK_EQUAL(spacer_pos, 8);
	}

BOOST_AUTO_TEST_SUITE_END()
