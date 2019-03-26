#include <mkl.h>

#include <iostream>
#include <vector>

int eig_dsyevd(double *matrix2eigvecs, double * eigvals, int L){
    //These nice values are inspired from armadillo. The prefactors give good performance.
    int lwork  = 2 * (1 + 6*L + 2*(L*L));
    int liwork = 3 * (3 + 5*L);
    int info   = 0;
    std::vector<double> work  ( lwork );
    std::vector<int   > iwork ( liwork );
    char jobz = 'V';
    info = LAPACKE_dsyevd_work(LAPACK_COL_MAJOR,jobz,'L',L,
                               matrix2eigvecs,
                               L,
                               eigvals,
                               work.data(),
                               lwork,
                               iwork.data(),
                               liwork);
    return info;
}


int main(){

    std::cout << "Working fine" << std::endl;

}