
#include <complex>


#include <Eigen/Core>
#include <h5pp/h5pp.h>
#include <math/class_eigsolver.h>


#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#ifdef MKL_AVAILABLE
#include <mkl_lapacke.h>
#include <mkl_service.h>
#include <mkl.h>
#else
#include <lapacke.h>
#endif



int eig_dsyevd(double *matrix2eigvecs, double * eigvals, int L){
    //These nice values are inspired from armadillo. The prefactors give good performance.
    int lwork  = 2 * (1 + 6*L + 2*(L*L));
    int liwork = 3 * (3 + 5*L);
    std::cout << " lwork  = " << lwork << std::endl;
    std::cout << " liwork = " << liwork << std::endl;

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

    return std::make_tuple(eigvals,eigvecs,info);
}


int main(){

    #ifdef _OPENMP
        omp_set_num_threads(std::thread::hardware_concurrency());
        Eigen::setNbThreads(std::thread::hardware_concurrency());
        std::cout << "Using Eigen  with " << Eigen::nbThreads() << " threads" << std::endl;
        std::cout << "Using OpenMP with " << omp_get_max_threads() << " threads" << std::endl;

        #ifdef OpenBLAS_AVAILABLE
        openblas_set_num_threads(std::thread::hardware_concurrency());
                    std::cout << OPENBLAS_VERSION
                              << " compiled with parallel mode " << openblas_get_parallel()
                              << " for target " << openblas_get_corename()
                              << " with config " << openblas_get_config()
                              << " with multithread threshold " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD
                              << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
        #endif

        #ifdef MKL_AVAILABLE
        mkl_set_num_threads(std::thread::hardware_concurrency());
        std::cout << "Using Intel MKL with " << mkl_get_max_threads() << " threads" << std::endl;
        #endif
    #endif



    Eigen::MatrixXd H_localA1, H_localA2;
    Eigen::MatrixXd H_localB1, H_localB2;

    Eigen::MatrixXd eigvecsA1, eigvecsA2;
    Eigen::MatrixXd eigvecsB1, eigvecsB2;
    Eigen::MatrixXd eigvalsA1, eigvalsA2;
    Eigen::MatrixXd eigvalsB1, eigvalsB2;
    Eigen::VectorXcd thetaA1, thetaA2;
    Eigen::VectorXcd thetaB1, thetaB2;

    Eigen::VectorXd overlapsA1, overlapsA2;
    Eigen::VectorXd overlapsB1, overlapsB2;

    h5pp::File fileA1(std::string(TEST_MATRIX_DIR) + "/testmatrices-A/lapacke_matrix-1.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileA2(std::string(TEST_MATRIX_DIR) + "/testmatrices-A/lapacke_matrix-2.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileB1(std::string(TEST_MATRIX_DIR) + "/testmatrices-B/lapacke_matrix-1.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    h5pp::File fileB2(std::string(TEST_MATRIX_DIR) + "/testmatrices-B/lapacke_matrix-2.h5", h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
    fileA1.readDataset(H_localA1,"matrix");
    fileA2.readDataset(H_localA2,"matrix");
    fileB1.readDataset(H_localB1,"matrix");
    fileB2.readDataset(H_localB2,"matrix");

    fileA1.readDataset(eigvecsA1,"eigvecs");
    fileA2.readDataset(eigvecsA2,"eigvecs");
    fileB1.readDataset(eigvecsB1,"eigvecs");
    fileB2.readDataset(eigvecsB2,"eigvecs");

    fileA1.readDataset(eigvalsA1,"eigvals");
    fileA2.readDataset(eigvalsA2,"eigvals");
    fileB1.readDataset(eigvalsB1,"eigvals");
    fileB2.readDataset(eigvalsB2,"eigvals");

    fileA1.readDataset(thetaA1,"theta");
    fileA2.readDataset(thetaA2,"theta");
    fileB1.readDataset(thetaB1,"theta");
    fileB2.readDataset(thetaB2,"theta");

    fileA1.readDataset(overlapsA1,"overlaps");
    fileA2.readDataset(overlapsA2,"overlaps");
    fileB1.readDataset(overlapsB1,"overlaps");
    fileB2.readDataset(overlapsB2,"overlaps");
    using namespace eigutils::eigSetting;

    class_eigsolver solverA;
    class_eigsolver solverB;

    [[maybe_unused]] auto [eigvalsA1_lapacke,eigvecsA1_lapacke,infoA1] = eig_dsyevd(H_localA1.data(),H_localA1.rows());
    [[maybe_unused]] auto [eigvalsB1_lapacke,eigvecsB1_lapacke,infoB1] = eig_dsyevd(H_localB1.data(),H_localB1.rows());

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


    auto overlap_A1A1_eig     =     (eigvecsA1.adjoint() * eigvecsA1_eig    ).diagonal().cwiseAbs().mean();
    auto overlap_A1B1_eig     =     (eigvecsA1.adjoint() * eigvecsB1_eig    ).diagonal().cwiseAbs().mean();
    auto overlap_A1A1_lapacke =     (eigvecsA1.adjoint() * eigvecsA1_lapacke).diagonal().cwiseAbs().mean();
    auto overlap_A1B1_lapacke =     (eigvecsA1.adjoint() * eigvecsB1_lapacke).diagonal().cwiseAbs().mean();

    [[maybe_unused]] auto [eigvalsA2_lapacke,eigvecsA2_lapacke,infoA2] = eig_dsyevd(H_localA2.data(),H_localA2.rows());
    [[maybe_unused]] auto [eigvalsB2_lapacke,eigvecsB2_lapacke,infoB2] = eig_dsyevd(H_localB2.data(),H_localB2.rows());
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

    auto overlap_A2A2_eig     =  (eigvecsA2.adjoint() * eigvecsA2_eig    ).diagonal().cwiseAbs().mean();
    auto overlap_A2B2_eig     =  (eigvecsA2.adjoint() * eigvecsB2_eig    ).diagonal().cwiseAbs().mean();
    auto overlap_A2A2_lapacke =  (eigvecsA2.adjoint() * eigvecsA2_lapacke).diagonal().cwiseAbs().mean();
    auto overlap_A2B2_lapacke =  (eigvecsA2.adjoint() * eigvecsB2_lapacke).diagonal().cwiseAbs().mean();




    std::cout << std::setprecision(16);
    std::cout << "H_localA1         vs H_localB1           diff   =   " << std::boolalpha << (H_localA1.array() - H_localB1.array()).cwiseAbs().sum() << std::endl;
    std::cout << "thetaA1           vs thetaB1             diff   =   " << std::boolalpha << (thetaA1.array()   - thetaB1.array()).cwiseAbs().sum() << std::endl;
    std::cout << "<thetaA1          |  thetaB1 >            diff  =   " << std::boolalpha << (thetaA1.adjoint() * thetaB1).cwiseAbs().sum() << std::endl;
    std::cout << "overlapsA1        vs overlapsB1          diff   =   " << std::boolalpha << (overlapsA1.array()- overlapsB1.array()).cwiseAbs().sum() << std::endl;
    std::cout << "eigvecsA1_read    vs eigvecsA1_lapacke   diff   =   " << std::boolalpha << (eigvecsA1.adjoint()         * eigvecsA1_lapacke).diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsB1_read    vs eigvecsB1_lapacke   diff   =   " << std::boolalpha << (eigvecsB1.adjoint()         * eigvecsB1_lapacke).diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsA1_read    vs eigvecsA1_solver    diff   =   " << std::boolalpha << (eigvecsA1.adjoint()         * eigvecsA1_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsB1_read    vs eigvecsB1_solver    diff   =   " << std::boolalpha << (eigvecsB1.adjoint()         * eigvecsB1_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsA1_lapacke vs eigvecsA1_solver    diff   =   " << std::boolalpha << (eigvecsA1_lapacke.adjoint() * eigvecsA1_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsB1_lapacke vs eigvecsB1_solver    diff   =   " << std::boolalpha << (eigvecsB1_lapacke.adjoint() * eigvecsB1_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvalsA1_read    vs eigvalsA1_lapacke   diff   =   " << std::boolalpha << (eigvalsA1.array()           - eigvalsA1_lapacke .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB1_read    vs eigvalsB1_lapacke   diff   =   " << std::boolalpha << (eigvalsB1.array()           - eigvalsB1_lapacke .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA1_read    vs eigvalsA1_solver    diff   =   " << std::boolalpha << (eigvalsA1.array()           - eigvalsA1_eig     .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB1_read    vs eigvalsB1_solver    diff   =   " << std::boolalpha << (eigvalsB1.array()           - eigvalsB1_eig     .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA1_lapacke vs eigvalsA1_solver    diff   =   " << std::boolalpha << (eigvalsA1_lapacke.array()   - eigvalsA1_eig     .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB1_lapacke vs eigvalsB1_solver    diff   =   " << std::boolalpha << (eigvalsB1_lapacke.array()   - eigvalsB1_eig     .array())  .cwiseAbs().sum() << std::endl;


    std::cout << "max_overlapA1_read     =  " << max_overlapA1_read     << std::endl;
    std::cout << "max_overlapA1_lapacke  =  " << max_overlapA1_lapacke  << std::endl;
    std::cout << "max_overlapA1_eig      =  " << max_overlapA1_eig      << std::endl;
    std::cout << "max_overlapB1_read     =  " << max_overlapB1_read     << std::endl;
    std::cout << "max_overlapB1_lapacke  =  " << max_overlapB1_lapacke  << std::endl;
    std::cout << "max_overlapB1_eig      =  " << max_overlapB1_eig      << std::endl;

    std::cout << "overlap_A1A1_eig       =  " << overlap_A1A1_eig       << std::endl;
    std::cout << "overlap_A1B1_eig       =  " << overlap_A1B1_eig       << std::endl;
    std::cout << "overlap_A1A1_lapacke   =  " << overlap_A1A1_lapacke   << std::endl;
    std::cout << "overlap_A1B1_lapacke   =  " << overlap_A1B1_lapacke   << std::endl;



    std::cout << "H_localA2         vs H_localB2           diff   =   " << std::boolalpha << (H_localA2.array() - H_localB2.array()).cwiseAbs().sum() << std::endl;
    std::cout << "thetaA2           vs thetaB2             diff   =   " << std::boolalpha << (thetaA2.array()   - thetaB2.array()).cwiseAbs().sum() << std::endl;
    std::cout << "<thetaA2          |  thetaB2 >            diff  =   " << std::boolalpha << (thetaA2.adjoint() * thetaB2).cwiseAbs().sum() << std::endl;
    std::cout << "overlapsA2        vs overlapsB2          diff   =   " << std::boolalpha << (overlapsA2.array()- overlapsB2.array()).cwiseAbs().sum() << std::endl;

    std::cout << "eigvecsA2_read    vs eigvecsA2_lapacke   diff   =   " << std::boolalpha << (eigvecsA2.adjoint()         * eigvecsA2_lapacke).diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsB2_read    vs eigvecsB2_lapacke   diff   =   " << std::boolalpha << (eigvecsB2.adjoint()         * eigvecsB2_lapacke).diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsA2_read    vs eigvecsA2_solver    diff   =   " << std::boolalpha << (eigvecsA2.adjoint()         * eigvecsA2_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsB2_read    vs eigvecsB2_solver    diff   =   " << std::boolalpha << (eigvecsB2.adjoint()         * eigvecsB2_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsA2_lapacke vs eigvecsA2_solver    diff   =   " << std::boolalpha << (eigvecsA2_lapacke.adjoint() * eigvecsA2_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvecsB2_lapacke vs eigvecsB2_solver    diff   =   " << std::boolalpha << (eigvecsB2_lapacke.adjoint() * eigvecsB2_eig)    .diagonal().cwiseAbs().mean() << std::endl;
    std::cout << "eigvalsA2_read    vs eigvalsA2_lapacke   diff   =   " << std::boolalpha << (eigvalsA2.array()           - eigvalsA2_lapacke .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB2_read    vs eigvalsB2_lapacke   diff   =   " << std::boolalpha << (eigvalsB2.array()           - eigvalsB2_lapacke .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA2_read    vs eigvalsA2_solver    diff   =   " << std::boolalpha << (eigvalsA2.array()           - eigvalsA2_eig     .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB2_read    vs eigvalsB2_solver    diff   =   " << std::boolalpha << (eigvalsB2.array()           - eigvalsB2_eig     .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsA2_lapacke vs eigvalsA2_solver    diff   =   " << std::boolalpha << (eigvalsA2_lapacke.array()   - eigvalsA2_eig     .array())  .cwiseAbs().sum() << std::endl;
    std::cout << "eigvalsB2_lapacke vs eigvalsB2_solver    diff   =   " << std::boolalpha << (eigvalsB2_lapacke.array()   - eigvalsB2_eig     .array())  .cwiseAbs().sum() << std::endl;


    std::cout << "max_overlapA2_read     =  " << max_overlapA2_read     << std::endl;
    std::cout << "max_overlapA2_lapacke  =  " << max_overlapA2_lapacke  << std::endl;
    std::cout << "max_overlapA2_eig      =  " << max_overlapA2_eig      << std::endl;
    std::cout << "max_overlapB2_read     =  " << max_overlapB2_read     << std::endl;
    std::cout << "max_overlapB2_lapacke  =  " << max_overlapB2_lapacke  << std::endl;
    std::cout << "max_overlapB2_eig      =  " << max_overlapB2_eig      << std::endl;

    std::cout << "overlap_A2A2_eig       =  " << overlap_A2A2_eig       << std::endl;
    std::cout << "overlap_A2B2_eig       =  " << overlap_A2B2_eig       << std::endl;
    std::cout << "overlap_A2A2_lapacke   =  " << overlap_A2A2_lapacke   << std::endl;
    std::cout << "overlap_A2B2_lapacke   =  " << overlap_A2B2_lapacke   << std::endl;

    if(std::abs(overlap_A1A1_eig     - 1) > 1e-10 ) throw std::runtime_error("overlap_A1A1_eig     : " + std::to_string(overlap_A1A1_eig     ));
    if(std::abs(overlap_A1B1_eig     - 1) > 1e-10 ) throw std::runtime_error("overlap_A1B1_eig     : " + std::to_string(overlap_A1B1_eig     ));
    if(std::abs(overlap_A1A1_lapacke - 1) > 1e-10 ) throw std::runtime_error("overlap_A1A1_lapacke : " + std::to_string(overlap_A1A1_lapacke ));
    if(std::abs(overlap_A1B1_lapacke - 1) > 1e-10 ) throw std::runtime_error("overlap_A1B1_lapacke : " + std::to_string(overlap_A1B1_lapacke ));



    if(std::abs(overlap_A2A2_eig     - 1) > 1e-10 ) throw std::runtime_error("overlap_A2A2_eig     : " + std::to_string(overlap_A2A2_eig     ));
    if(std::abs(overlap_A2B2_eig     - 1) > 1e-10 ) throw std::runtime_error("overlap_A2B2_eig     : " + std::to_string(overlap_A2B2_eig     ));
    if(std::abs(overlap_A2A2_lapacke - 1) > 1e-10 ) throw std::runtime_error("overlap_A2A2_lapacke : " + std::to_string(overlap_A2A2_lapacke ));
    if(std::abs(overlap_A2B2_lapacke - 1) > 1e-10 ) throw std::runtime_error("overlap_A2B2_lapacke : " + std::to_string(overlap_A2B2_lapacke ));


//    std::cout << "theta * eigvecsA2         = \n" << std::fixed << (thetaA2.adjoint() * eigvecsA2).cwiseAbs().transpose() << std::endl;


}