/*
    Fixes a compiler error caused by a duplicate symbol in lapack and boost.
    See: https://github.com/xianyi/OpenBLAS/issues/1992
*/

#include <complex>
#undef lapack_complex_float
#undef lapack_complex_double
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
