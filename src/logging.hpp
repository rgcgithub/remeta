#ifndef REMETA_LOGGING_H
#define REMETA_LOGGING_H

// #define BOOST_LOG_DYN_LINK 1

#include <string>
using namespace std;

void init_logging(string loglevel, string logfile);
void log_debug(const string& s);
void log_info(const string& s);
void log_warning(const string& s);
void log_error(const string& s, const int& die = 0);

#endif
