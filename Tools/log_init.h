#pragma once

#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>

namespace logging = boost::log;
namespace src = logging::sources;
namespace keywords = logging::keywords;

static void init_log(bool verbose, bool debug, const std::string &main_file_name = "main.log",
					 const std::string &debug_file_name = "debug.log")
{
	typedef boost::shared_ptr<logging::sinks::synchronous_sink<logging::sinks::text_ostream_backend> > console_sink_t;
	typedef boost::shared_ptr<logging::sinks::synchronous_sink<logging::sinks::text_file_backend> > file_sink_t;

	logging::trivial::severity_level min_level = verbose ? logging::trivial::trace
														 : (debug ? logging::trivial::debug
																  : logging::trivial::info);
	logging::core::get()->set_filter(
			logging::trivial::severity >= min_level
	);

	logging::formatter logFmt =
			logging::expressions::format("%1%")
			% logging::expressions::smessage;


	console_sink_t console_sink = logging::add_console_log(std::cerr);
	console_sink->set_filter(
			logging::trivial::severity >= logging::trivial::info ||
			logging::trivial::severity == logging::trivial::trace
	);
	console_sink->set_formatter(logFmt);

	file_sink_t fs_sink = logging::add_file_log(
			logging::keywords::file_name = main_file_name,
			logging::keywords::open_mode = std::ios_base::out,
			logging::keywords::filter = (logging::trivial::severity >= logging::trivial::info ||
										 logging::trivial::severity == logging::trivial::trace));
	fs_sink->set_formatter(logFmt);
	fs_sink->locked_backend()->auto_flush(true);

	if (debug)
	{
		file_sink_t debug_fs_sink = logging::add_file_log(
				logging::keywords::file_name = debug_file_name,
				logging::keywords::open_mode = std::ios_base::out,
				logging::keywords::filter = (logging::trivial::severity >= logging::trivial::debug));
		debug_fs_sink->set_formatter(logFmt);
		debug_fs_sink->locked_backend()->auto_flush(true);
	}
}

static void init_test_logs()
{
	logging::core::get()->set_filter(
			logging::trivial::severity >= logging::trivial::fatal
	);
}