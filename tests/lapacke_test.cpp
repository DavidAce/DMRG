#include <Eigen/Core>
#include <mps_routines/nmspc_mps_tools.h>
#include <complex>
#include <mps_routines/class_superblock.h>
#include <LBFGS.h>
#include <h5pp/h5pp.h>
#include <general/class_eigsolver.h>
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

    Eigen::VectorXcd thetaA0,thetaA1, thetaA2;
    Eigen::VectorXcd thetaB0,thetaB1, thetaB2;

    Eigen::VectorXd overlapsA0,overlapsA1, overlapsA2;
    Eigen::VectorXd overlapsB0,overlapsB1, overlapsB2;

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

    fileA0.readDataset(eigvecsA0,"eigvecs");
    fileA1.readDataset(eigvecsA1,"eigvecs");
    fileA2.readDataset(eigvecsA2,"eigvecs");
    fileB0.readDataset(eigvecsB0,"eigvecs");
    fileB1.readDataset(eigvecsB1,"eigvecs");
    fileB2.readDataset(eigvecsB2,"eigvecs");

    fileA0.readDataset(eigvalsA0,"eigvals");
    fileA1.readDataset(eigvalsA1,"eigvals");
    fileA2.readDataset(eigvalsA2,"eigvals");
    fileB0.readDataset(eigvalsB0,"eigvals");
    fileB1.readDataset(eigvalsB1,"eigvals");
    fileB2.readDataset(eigvalsB2,"eigvals");

    fileA0.readDataset(thetaA0,"theta");
    fileA1.readDataset(thetaA1,"theta");
    fileA2.readDataset(thetaA2,"theta");
    fileB0.readDataset(thetaB0,"theta");
    fileB1.readDataset(thetaB1,"theta");
    fileB2.readDataset(thetaB2,"theta");

    fileA0.readDataset(overlapsA0,"overlaps");
    fileA1.readDataset(overlapsA1,"overlaps");
    fileA2.readDataset(overlapsA2,"overlaps");
    fileB0.readDataset(overlapsB0,"overlaps");
    fileB1.readDataset(overlapsB1,"overlaps");
    fileB2.readDataset(overlapsB2,"overlaps");
    using namespace eigutils::eigSetting;

    class_eigsolver solverA;
    class_eigsolver solverB;


    auto [eigvalsA0_lapacke,eigvecsA0_lapacke,infoA0] = eig_dsyevd(H_localA0.data(),H_localA0.rows());
    auto [eigvalsB0_lapacke,eigvecsB0_lapacke,infoB0] = eig_dsyevd(H_localB0.data(),H_localB0.rows());

    solverA.eig<Type::REAL, Form::SYMMETRIC>(H_localA0,true,false);
    solverB.eig<Type::REAL, Form::SYMMETRIC>(H_localB0,true,false);
    Eigen::VectorXd eigvalsA0_eig           = Eigen::Map<const Eigen::VectorXd> (solverA.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solverA.solution.meta.cols);
    Eigen::MatrixXd eigvecsA0_eig           = Eigen::Map<const Eigen::MatrixXd> (solverA.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solverA.solution.meta.rows,solverA.solution.meta.cols);
    Eigen::VectorXd eigvalsB0_eig           = Eigen::Map<const Eigen::VectorXd> (solverB.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solverB.solution.meta.cols);
    Eigen::MatrixXd eigvecsB0_eig           = Eigen::Map<const Eigen::MatrixXd> (solverB.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solverB.solution.meta.rows,solverB.solution.meta.cols);

    Eigen::VectorXd overlapsA0_read         = (thetaA0.adjoint() * eigvecsA0).cwiseAbs();
    Eigen::VectorXd overlapsA0_lapacke      = (thetaA0.adjoint() * eigvecsA0_lapacke).cwiseAbs();
    Eigen::VectorXd overlapsA0_eig          = (thetaA0.adjoint() * eigvecsA0_eig).cwiseAbs();
    Eigen::VectorXd overlapsB0_read         = (thetaB0.adjoint() * eigvecsB0).cwiseAbs();
    Eigen::VectorXd overlapsB0_lapacke      = (thetaB0.adjoint() * eigvecsB0_lapacke).cwiseAbs();
    Eigen::VectorXd overlapsB0_eig          = (thetaB0.adjoint() * eigvecsB0_eig).cwiseAbs();


    double max_overlapA0_read        = overlapsA0_read    .maxCoeff();
    double max_overlapA0_lapacke     = overlapsA0_lapacke .maxCoeff();
    double max_overlapA0_eig         = overlapsA0_eig     .maxCoeff();
    double max_overlapB0_read        = overlapsB0_read    .maxCoeff();
    double max_overlapB0_lapacke     = overlapsB0_lapacke .maxCoeff();
    double max_overlapB0_eig         = overlapsB0_eig     .maxCoeff();






    auto [eigvalsA1_lapacke,eigvecsA1_lapacke,infoA1] = eig_dsyevd(H_localA1.data(),H_localA1.rows());
    auto [eigvalsB1_lapacke,eigvecsB1_lapacke,infoB1] = eig_dsyevd(H_localB1.data(),H_localB1.rows());

    solverA.eig<Type::REAL, Form::SYMMETRIC>(H_localA1,true,false);
    solverB.eig<Type::REAL, Form::SYMMETRIC>(H_localB1,true,false);
    Eigen::VectorXd eigvalsA1_eig           = Eigen::Map<const Eigen::VectorXd> (solverA.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solverA.solution.meta.cols);
    Eigen::MatrixXd eigvecsA1_eig           = Eigen::Map<const Eigen::MatrixXd> (solverA.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solverA.solution.meta.rows,solverA.solution.meta.cols);
    Eigen::VectorXd eigvalsB1_eig           = Eigen::Map<const Eigen::VectorXd> (solverB.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solverB.solution.meta.cols);
    Eigen::MatrixXd eigvecsB1_eig           = Eigen::Map<const Eigen::MatrixXd> (solverB.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solverB.solution.meta.rows,solverB.solution.meta.cols);

    Eigen::VectorXd overlapsA1_read         = (thetaA1.adjoint() * eigvecsA1).cwiseAbs();
    Eigen::VectorXd overlapsA1_lapacke      = (thetaA1.adjoint() * eigvecsA1_lapacke).cwiseAbs();
    Eigen::VectorXd overlapsA1_eig          = (thetaA1.adjoint() * eigvecsA1_eig).cwiseAbs();
    Eigen::VectorXd overlapsB1_read         = (thetaB1.adjoint() * eigvecsB1).cwiseAbs();
    Eigen::VectorXd overlapsB1_lapacke      = (thetaB1.adjoint() * eigvecsB1_lapacke).cwiseAbs();
    Eigen::VectorXd overlapsB1_eig          = (thetaB1.adjoint() * eigvecsB1_eig).cwiseAbs();


    double max_overlapA1_read        = overlapsA1_read    .maxCoeff();
    double max_overlapA1_lapacke     = overlapsA1_lapacke .maxCoeff();
    double max_overlapA1_eig         = overlapsA1_eig     .maxCoeff();
    double max_overlapB1_read        = overlapsB1_read    .maxCoeff();
    double max_overlapB1_lapacke     = overlapsB1_lapacke .maxCoeff();
    double max_overlapB1_eig         = overlapsB1_eig     .maxCoeff();







    auto [eigvalsA2_lapacke,eigvecsA2_lapacke,infoA2] = eig_dsyevd(H_localA2.data(),H_localA2.rows());
    auto [eigvalsB2_lapacke,eigvecsB2_lapacke,infoB2] = eig_dsyevd(H_localB2.data(),H_localB2.rows());
    solverA.eig<Type::REAL, Form::SYMMETRIC>(H_localA2,true,false);
    solverB.eig<Type::REAL, Form::SYMMETRIC>(H_localB2,true,false);
    Eigen::VectorXd eigvalsA2_eig           = Eigen::Map<const Eigen::VectorXd> (solverA.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solverA.solution.meta.cols);
    Eigen::MatrixXd eigvecsA2_eig           = Eigen::Map<const Eigen::MatrixXd> (solverA.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solverA.solution.meta.rows,solverA.solution.meta.cols);
    Eigen::VectorXd eigvalsB2_eig           = Eigen::Map<const Eigen::VectorXd> (solverB.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solverB.solution.meta.cols);
    Eigen::MatrixXd eigvecsB2_eig           = Eigen::Map<const Eigen::MatrixXd> (solverB.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solverB.solution.meta.rows,solverB.solution.meta.cols);

    Eigen::VectorXd overlapsA2_read         = (thetaA2.adjoint() * eigvecsA2).cwiseAbs();
    Eigen::VectorXd overlapsA2_lapacke      = (thetaA2.adjoint() * eigvecsA2_lapacke).cwiseAbs();
    Eigen::VectorXd overlapsA2_eig          = (thetaA2.adjoint() * eigvecsA2_eig).cwiseAbs();
    Eigen::VectorXd overlapsB2_read         = (thetaB2.adjoint() * eigvecsB2).cwiseAbs();
    Eigen::VectorXd overlapsB2_lapacke      = (thetaB2.adjoint() * eigvecsB2_lapacke).cwiseAbs();
    Eigen::VectorXd overlapsB2_eig          = (thetaB2.adjoint() * eigvecsB2_eig).cwiseAbs();


    double max_overlapA2_read        = overlapsA2_read    .maxCoeff();
    double max_overlapA2_lapacke     = overlapsA2_lapacke .maxCoeff();
    double max_overlapA2_eig         = overlapsA2_eig     .maxCoeff();
    double max_overlapB2_read        = overlapsB2_read    .maxCoeff();
    double max_overlapB2_lapacke     = overlapsB2_lapacke .maxCoeff();
    double max_overlapB2_eig         = overlapsB2_eig     .maxCoeff();


    std::cout << std::setprecision(16);

    std::cout << "H_localA0         vs H_localB0           diff   =   " << std::boolalpha << (H_localA0.array() - H_localB0.array()).cwiseAbs().sum() << std::endl;
    std::cout << "thetaA0           vs thetaB0             diff   =   " << std::boolalpha << (thetaA0.array()   - thetaB0.array()).cwiseAbs().sum() << std::endl;
    std::cout << "overlapsA0        vs overlapsB0          diff   =   " << std::boolalpha << (overlapsA0.array()- overlapsB0.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA0_read    vs eigvalsA0_lapacke   diff   =   " << std::boolalpha << (eigvalsA0.array() - eigvalsA0_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA0_read    vs eigvecsA0_lapacke   diff   =   " << std::boolalpha << (eigvecsA0.array() - eigvecsA0_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB0_read    vs eigvalsB0_lapacke   diff   =   " << std::boolalpha << (eigvalsB0.array() - eigvalsB0_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB0_read    vs eigvecsB0_lapacke   diff   =   " << std::boolalpha << (eigvecsB0.array() - eigvecsB0_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA0_read    vs eigvalsA0_solver    diff   =   " << std::boolalpha << (eigvalsA0.array() - eigvalsA0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA0_read    vs eigvecsA0_solver    diff   =   " << std::boolalpha << (eigvecsA0.array() - eigvecsA0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB0_read    vs eigvalsB0_solver    diff   =   " << std::boolalpha << (eigvalsB0.array() - eigvalsB0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB0_read    vs eigvecsB0_solver    diff   =   " << std::boolalpha << (eigvecsB0.array() - eigvecsB0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA0_lapacke vs eigvalsA0_solver    diff   =   " << std::boolalpha << (eigvalsA0_lapacke.array() - eigvalsA0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA0_lapacke vs eigvecsA0_solver    diff   =   " << std::boolalpha << (eigvecsA0_lapacke.array() - eigvecsA0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB0_lapacke vs eigvalsB0_solver    diff   =   " << std::boolalpha << (eigvalsB0_lapacke.array() - eigvalsB0_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB0_lapacke vs eigvecsB0_solver    diff   =   " << std::boolalpha << (eigvecsB0_lapacke.array() - eigvecsB0_eig.array()).cwiseAbs().sum() << std::endl;

    std::cout << "max_overlapA0_read     =  " << max_overlapA0_read     << std::endl;
    std::cout << "max_overlapA0_lapacke  =  " << max_overlapA0_lapacke  << std::endl;
    std::cout << "max_overlapA0_eig      =  " << max_overlapA0_eig      << std::endl;
    std::cout << "max_overlapB0_read     =  " << max_overlapB0_read     << std::endl;
    std::cout << "max_overlapB0_lapacke  =  " << max_overlapB0_lapacke  << std::endl;
    std::cout << "max_overlapB0_eig      =  " << max_overlapB0_eig      << std::endl;



    std::cout << "H_localA1         vs H_localB1           diff   =   " << std::boolalpha << (H_localA1.array() - H_localB1.array()).cwiseAbs().sum() << std::endl;
    std::cout << "thetaA1           vs thetaB1             diff   =   " << std::boolalpha << (thetaA1.array()   - thetaB1.array()).cwiseAbs().sum() << std::endl;
    std::cout << "overlapsA1        vs overlapsB1          diff   =   " << std::boolalpha << (overlapsA1.array()- overlapsB1.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA1_read    vs eigvalsA1_lapacke   diff   =   " << std::boolalpha << (eigvalsA1.array() - eigvalsA1_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA1_read    vs eigvecsA1_lapacke   diff   =   " << std::boolalpha << (eigvecsA1.array() - eigvecsA1_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB1_read    vs eigvalsB1_lapacke   diff   =   " << std::boolalpha << (eigvalsB1.array() - eigvalsB1_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB1_read    vs eigvecsB1_lapacke   diff   =   " << std::boolalpha << (eigvecsB1.array() - eigvecsB1_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA1_read    vs eigvalsA1_solver    diff   =   " << std::boolalpha << (eigvalsA1.array() - eigvalsA1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA1_read    vs eigvecsA1_solver    diff   =   " << std::boolalpha << (eigvecsA1.array() - eigvecsA1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB1_read    vs eigvalsB1_solver    diff   =   " << std::boolalpha << (eigvalsB1.array() - eigvalsB1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB1_read    vs eigvecsB1_solver    diff   =   " << std::boolalpha << (eigvecsB1.array() - eigvecsB1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA1_lapacke vs eigvalsA1_solver    diff   =   " << std::boolalpha << (eigvalsA1_lapacke.array() - eigvalsA1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA1_lapacke vs eigvecsA1_solver    diff   =   " << std::boolalpha << (eigvecsA1_lapacke.array() - eigvecsA1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB1_lapacke vs eigvalsB1_solver    diff   =   " << std::boolalpha << (eigvalsB1_lapacke.array() - eigvalsB1_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB1_lapacke vs eigvecsB1_solver    diff   =   " << std::boolalpha << (eigvecsB1_lapacke.array() - eigvecsB1_eig.array()).cwiseAbs().sum() << std::endl;

    std::cout << "max_overlapA1_read     =  " << max_overlapA1_read     << std::endl;
    std::cout << "max_overlapA1_lapacke  =  " << max_overlapA1_lapacke  << std::endl;
    std::cout << "max_overlapA1_eig      =  " << max_overlapA1_eig      << std::endl;
    std::cout << "max_overlapB1_read     =  " << max_overlapB1_read     << std::endl;
    std::cout << "max_overlapB1_lapacke  =  " << max_overlapB1_lapacke  << std::endl;
    std::cout << "max_overlapB1_eig      =  " << max_overlapB1_eig      << std::endl;



    std::cout << "H_localA2         vs H_localB2           diff   =   " << std::boolalpha << (H_localA2.array() - H_localB2.array()).cwiseAbs().sum() << std::endl;
    std::cout << "thetaA2           vs thetaB2             diff   =   " << std::boolalpha << (thetaA2.array()   - thetaB2.array()).cwiseAbs().sum() << std::endl;
    std::cout << "overlapsA2        vs overlapsB2          diff   =   " << std::boolalpha << (overlapsA2.array()- overlapsB2.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA2_read    vs eigvalsA2_lapacke   diff   =   " << std::boolalpha << (eigvalsA2.array() - eigvalsA2_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA2_read    vs eigvecsA2_lapacke   diff   =   " << std::boolalpha << (eigvecsA2.array() - eigvecsA2_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB2_read    vs eigvalsB2_lapacke   diff   =   " << std::boolalpha << (eigvalsB2.array() - eigvalsB2_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB2_read    vs eigvecsB2_lapacke   diff   =   " << std::boolalpha << (eigvecsB2.array() - eigvecsB2_lapacke.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA2_read    vs eigvalsA2_solver    diff   =   " << std::boolalpha << (eigvalsA2.array() - eigvalsA2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA2_read    vs eigvecsA2_solver    diff   =   " << std::boolalpha << (eigvecsA2.array() - eigvecsA2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB2_read    vs eigvalsB2_solver    diff   =   " << std::boolalpha << (eigvalsB2.array() - eigvalsB2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB2_read    vs eigvecsB2_solver    diff   =   " << std::boolalpha << (eigvecsB2.array() - eigvecsB2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA2_lapacke vs eigvalsA2_solver    diff   =   " << std::boolalpha << (eigvalsA2_lapacke.array() - eigvalsA2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA2_lapacke vs eigvecsA2_solver    diff   =   " << std::boolalpha << (eigvecsA2_lapacke.array() - eigvecsA2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB2_lapacke vs eigvalsB2_solver    diff   =   " << std::boolalpha << (eigvalsB2_lapacke.array() - eigvalsB2_eig.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsB2_lapacke vs eigvecsB2_solver    diff   =   " << std::boolalpha << (eigvecsB2_lapacke.array() - eigvecsB2_eig.array()).cwiseAbs().sum() << std::endl;

    std::cout << "max_overlapA2_read     =  " << max_overlapA2_read     << std::endl;
    std::cout << "max_overlapA2_lapacke  =  " << max_overlapA2_lapacke  << std::endl;
    std::cout << "max_overlapA2_eig      =  " << max_overlapA2_eig      << std::endl;
    std::cout << "max_overlapB2_read     =  " << max_overlapB2_read     << std::endl;
    std::cout << "max_overlapB2_lapacke  =  " << max_overlapB2_lapacke  << std::endl;
    std::cout << "max_overlapB2_eig      =  " << max_overlapB2_eig      << std::endl;

//    std::cout << "theta * eigvecsA2         = \n" << std::fixed << (thetaA2.adjoint() * eigvecsA2).cwiseAbs().transpose() << std::endl;


}