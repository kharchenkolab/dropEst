#pragma once

#include <boost/log/trivial.hpp>

#define L_DEBUG BOOST_LOG_TRIVIAL(debug)
#define L_INFO BOOST_LOG_TRIVIAL(info)
#define L_ERR BOOST_LOG_TRIVIAL(error)
