#include "logging.hpp"

#include <iostream>
#include <string>
using namespace std;

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/support/date_time.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;

void init_logging(string loglevel, string logfile) {
  if (logfile != "") {
    logging::add_file_log(
      keywords::file_name = logfile,
      keywords::auto_flush = true
    );
  }
  logging::add_console_log(
    std::cout,
    keywords::filter = logging::trivial::severity < logging::trivial::warning
  );
  logging::add_console_log(
    std::cerr,
    keywords::filter = logging::trivial::severity >= logging::trivial::warning
  );

  if (loglevel == "debug") {
    logging::core::get()->set_filter(
      logging::trivial::severity >= logging::trivial::debug
    );
  } else {
    logging::core::get()->set_filter(
      logging::trivial::severity >= logging::trivial::info
    );
  }
}

void log_debug(const string& s) {
  BOOST_LOG_TRIVIAL(debug) << "debug: " << s;
}

void log_info(const string& s) {
  BOOST_LOG_TRIVIAL(info) << s;
}

void log_warning(const string& s) {
  BOOST_LOG_TRIVIAL(warning) << "warning: " << s;
}

void log_error(const string& s, const int& die) {
  BOOST_LOG_TRIVIAL(error) << "error: " << s;
  if (die) {
    exit(die);
  }
}