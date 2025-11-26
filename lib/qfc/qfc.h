// Code was obtained from Robert Davies' webpage (http://www.robertnz.net/download.html) and was modified to remove the use of functions from setjmp library which does not work well with C++.

#ifndef QFC_H
#define QFC_H

//#define UseDouble 0             /* all floating point double */

#ifdef __cplusplus
extern "C" {
#endif

double qf(double*, double*, int*, int, double, double, int, double, double*,
          int*);

#ifdef __cplusplus
}
#endif

#endif
