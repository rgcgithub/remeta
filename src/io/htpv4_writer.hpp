/* htpv4_writer.hpp
* Author: Tyler Joseph
* 
* Write raw or bgzipped HTPv4 files.
* 
* Inherits from BgzWriter class.
* 
* Example:
*   // write a bgzipped file
*   HTPv4Writer writer("myfile.gz", "w");
*   writer.write_header();
*   writer.writerec(htpv4_record_t);
*   writer.write(string);
*
*   // write an uncompressed file file
*   HTPv4Writer writer("myfile", "wu");
* 
*   // append to a file
*   HTPv4Writer writer("filewithdata.gz", "w");
* 
*   // close the file
*   writer.close();
*/
#ifndef HTPv4_WRITER_H
#define HTPv4_WRITER_H

#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
using namespace std;

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include "bgz_writer.hpp"
#include "htpv4_reader.hpp"

class HTPv4Writer : public BgzWriter {
  public:
    HTPv4Writer() : BgzWriter(){};
    HTPv4Writer(string filepath, string mode) : BgzWriter(filepath, mode){};
    ~HTPv4Writer(){};
    void writerec(htpv4_record_t rec);
    void writeheader();
};

#endif