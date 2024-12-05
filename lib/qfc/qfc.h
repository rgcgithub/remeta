// Code was obtained from Robert Davies' webpage (http://www.robertnz.net/download.html) and was modified to remove the use of functions from setjmp library which does not work well with C++.

#ifndef QFC_H
#define QFC_H

//#define UseDouble 0             /* all floating point double */

#define TRUE 1
#define FALSE 0
typedef int BOOL;

#define pi 3.14159265358979
#define log28 .0866 /*  log(2.0) / 8.0  */

#ifdef __cplusplus
extern "C" {
#endif

double qf(double*, double*, int*, int, double, double, int, double, double*,
          int*);

#ifdef __cplusplus
}
#endif

#endif
