#pragma once

#include <string>

#include <boost/log/trivial.hpp>

#define L_TRACE BOOST_LOG_TRIVIAL(trace)
#define L_DEBUG BOOST_LOG_TRIVIAL(debug)
#define L_INFO BOOST_LOG_TRIVIAL(info)
#define L_ERR BOOST_LOG_TRIVIAL(error)

namespace Tools
{
	void init_log(bool verbose, bool debug, const std::string &main_file_name = "main.log",
						 const std::string &debug_file_name = "debug.log");

	void init_test_logs(boost::log::trivial::severity_level level = boost::log::trivial::fatal);
}