#include <complex>
#ifndef lapack_complex_float
    #define lapack_complex_float  std::complex<float>
#endif
#ifndef lapack_complex_double
    #define lapack_complex_double std::complex<double>
#endif
#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#elif __has_include(<openblas/lapacke.h>)
#include <openblas/lapacke.h>
#else
#include <lapacke.h>
#endif

int main (int argc, const char * argv[]){
    double a[5][3] = {{1,1,1},{2,3,4},{3,5,2},{4,2,5},{5,4,3}};
    double b[5][2] = {{-10,-3},{12,14},{14,12},{16,16},{18,16}};
    lapack_int info,m,n,lda,ldb,nrhs;
    int i,j;
    m = 5;
    n = 3;
    nrhs = 2;
    lda = 3;
    ldb = 2;
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
    return(info);
}
