#pragma once

#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>

namespace logging = boost::log;
namespace src = logging::sources;
namespace keywords = logging::keywords;

static void init_log(logging::trivial::severity_level min_level)
{
	logging::core::get()->set_filter(
			logging::trivial::severity >= min_level
	);

	logging::formatter logFmt =
			logging::expressions::format("%1%")
			% logging::expressions::smessage;

	auto consoleSink = logging::add_console_log(std::cerr, keywords::filter = (logging::trivial::severity >=
																			   logging::trivial::info));
	consoleSink->set_formatter(logFmt);

	auto fsSink = logging::add_file_log(
			logging::keywords::file_name = "debug.log",
			logging::keywords::open_mode = std::ios_base::out,
			logging::keywords::filter = (logging::trivial::severity >= logging::trivial::debug));
	fsSink->set_formatter(logFmt);
	fsSink->locked_backend()->auto_flush(true);
}