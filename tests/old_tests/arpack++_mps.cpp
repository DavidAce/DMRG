#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#include <openblas_config.h>
#endif

#ifdef MKL_AVAILABLE
#include <mkl_service.h>
#include <mkl.h>
#endif




#include <math/class_eigsolver.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>


using namespace eigutils::eigSetting;

// Find the optimal MPS representation given an effective Hamiltonian made of left, right environments and a Hamiltonian MPO
// called  L,R and M respectively. Start from an initial guess theta.


int main()
{

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

    using Scalar = std::complex<double>;
    int nev = 1;
    int eigMaxNcv = 16;
    std::array<long,4> shape_theta4 = {2,1,2,1};
    std::array<long,4> shape_mpo4   = {3,3,2,2};
    std::vector<Scalar> L = {0.924445, 0.381316, 1.0};
    std::vector<Scalar> R = {1.0,  0.356012 , -0.934481};
    std::vector<Scalar> M = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             -1.0, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<Scalar> theta = {-0.191154, 0.964728,
                                 -0.0351791,0.177544};

    DenseHamiltonianProduct<Scalar> matrix(L.data(),R.data(),M.data(),M.data(),shape_theta4,shape_mpo4);
    class_eigsolver solver;
    solver.eigs_dense(matrix,nev,eigMaxNcv, NAN,Form::SYMMETRIC,Ritz::LM,Side::R, true,false);

    using namespace Textra;
    using namespace Eigen;
    [[maybe_unused]] auto eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solver.solution.meta.cols);
    [[maybe_unused]] auto eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
//
//    Eigen::TensorMap<const Eigen::Tensor<const Scalar,2>> eigvecs (solver.ref_eigvecs().data(), shape_theta4[0]*shape_theta4[1], shape_theta4[2]*shape_theta4[3]);
//    Eigen::TensorMap<const Eigen::Tensor<const Scalar,1>> eigvals (solver.ref_eigvals().data(), nev);


    double E = std::real(eigvals(0))/2.0;
    std::cout << "Energy: " << std::setprecision(10) << E << std::endl;
    std::cout << "Diff  : " << std::abs(E - -1.110756115) << std::endl;
    if (std::abs(E - -1.110756115) > 1e-6){return 1;}
    else {return 0;}
}