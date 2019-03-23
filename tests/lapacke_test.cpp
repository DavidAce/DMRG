#include <Eigen/Core>

//#include <general/class_eigsolver.h>
#include <h5pp/h5pp.h>
#include <lapacke.h>


int eig_dsyevd(double *matrix2eigvecs, double * eigvals, int L){
    //These nice values are inspired from armadillo. The prefactors give good performance.
    int lwork  =  2 * (1 + 6*L + 2*(L*L));
    int liwork =  3 * (3 + 5*L);
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

std::tuple<Eigen::VectorXd,Eigen::MatrixXd, int> eig_dsyevd(const double* matrix, int L){
//    using namespace eigutils::eigSetting;
    Eigen::VectorXd eigvals(L);
    Eigen::MatrixXd eigvecs(L,L);
    std::copy(matrix, matrix + L*L, eigvecs.data());
    int info = eig_dsyevd(eigvecs.data(),eigvals.data(), L);
//    if (info == 0){
//        solution.meta.eigvecs_found = true;
//        solution.meta.eigvals_found = true;
//        solution.meta.rows           = L;
//        solution.meta.cols           = L;
//        solution.meta.nev            = L;
//        solution.meta.n              = L;
//        solution.meta.form           = Form::SYMMETRIC;
//        solution.meta.type           = Type::REAL ;
//        int i = 0;
//        int j = 0;
//        int c = 0;
//        if(L == 64){
//            for(auto & val : eigvecs){
//                std::cout << "eigvec(" << i << "," << j << ") = " << std::setprecision(16) << std::setw(20) <<  std::fixed<< std::right  << val
//                          <<  " matrix = "  << std::setw(20) << std::fixed << std::right  << matrix[c] << std::endl;
//                i++;
//                c++;
//                if (i == L){i = 0; j++;}
//            }
//        }

//    }else{
//        throw std::runtime_error("LAPACK dsyevd failed with error: " + std::to_string(info));
//    }


    return std::make_tuple(eigvals,eigvecs,info);
}


int main(){

    Eigen::MatrixXd H_localA0,H_localA1, H_localA2;
    Eigen::MatrixXd H_localB0,H_localB1, H_localB2;

    Eigen::MatrixXd eigvecsA0,eigvecsA1, eigvecsA2;
    Eigen::MatrixXd eigvecsB0,eigvecsB1, eigvecsB2;

    Eigen::MatrixXd eigvalsA0,eigvalsA1, eigvalsA2;
    Eigen::MatrixXd eigvalsB0,eigvalsB1, eigvalsB2;

    h5pp::File fileA0("../tests/testmatrices-A/lapacke_matrix.h5"  , h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileA1("../tests/testmatrices-A/lapacke_matrix-1.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileA2("../tests/testmatrices-A/lapacke_matrix-2.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileB0("../tests/testmatrices-B/lapacke_matrix.h5"  , h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileB1("../tests/testmatrices-B/lapacke_matrix-1.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileB2("../tests/testmatrices-B/lapacke_matrix-2.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    fileA0.readDataset(H_localA0,"matrix");
    fileA1.readDataset(H_localA1,"matrix");
    fileA2.readDataset(H_localA2,"matrix");
    fileB0.readDataset(H_localB0,"matrix");
    fileB1.readDataset(H_localB1,"matrix");
    fileB2.readDataset(H_localB2,"matrix");

    std::cout << std::boolalpha << (H_localA2.array() - H_localB2.array()).cwiseAbs().sum() << std::endl;

    auto [eigvalsA2_calc,eigvecsA2_calc,infoA2] = eig_dsyevd(H_localA2.data(),H_localA2.rows());
    auto [eigvalsB2_calc,eigvecsB2_calc,infoB2] = eig_dsyevd(H_localB2.data(),H_localB2.rows());
    std::cout << std::boolalpha << (eigvalsA2_calc.array() - eigvalsB2_calc.array()).cwiseAbs().sum() << std::endl;
    std::cout << std::boolalpha << (eigvecsA2_calc.array() - eigvecsB2_calc.array()).cwiseAbs().sum() << std::endl;




}