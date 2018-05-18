
#include <general/class_arpack_eigsolver.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
using namespace std::complex_literals;
// Find the optimal MPS representation given an effective Hamiltonian made of left, right environments and a Hamiltonian MPO
// called  L,R and M respectively. Start from an initial guess theta.


int main()
{
    using Scalar = std::complex<double>;
    int nev = 1;
    int ncv = 16;
    int maxiter = 1000;
    double threshold = 1e-12;
    std::array<long,4> shape_theta4 = {2,1,2,1};
    std::array<long,4> shape_mpo4   = {3,3,2,2};
    std::vector<Scalar> L = {0.924445, 0.381316, 1.0};
    std::vector<Scalar> R = {1.0,  0.356012 , -0.934481};
    std::vector<std::complex<double>> M = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                           -1.0, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<Scalar> theta = {-0.191154, 0.964728,
                                 -0.0351791,0.177544};

    class_arpack_eigsolver<Scalar, Form::GENERAL> solver(L.data(),
                                                         R.data(),
                                                         M.data(),
                                                         M.data(),
                                                         shape_theta4,
                                                         shape_mpo4,
                                                         Ritz::SR,
                                                         nev,
                                                         ncv,
                                                         threshold,
                                                         maxiter,
                                                         true ,
                                                         theta.data());
    using namespace Textra;
    using namespace Eigen;
    Eigen::TensorMap<const Eigen::Tensor<const Scalar,2>> eigvecs (solver.ref_eigvecs().data(), shape_theta4[0]*shape_theta4[1], shape_theta4[2]*shape_theta4[3]);
    Eigen::TensorMap<const Eigen::Tensor<const Scalar,1>> eigvals (solver.ref_eigvals().data(), nev);


    double E = std::real(eigvals(0))/2.0;
    std::cout << "Energy: " << std::setprecision(10) << E << std::endl;
    std::cout << "Diff  : " << std::abs(E - -1.110756115) << std::endl;
    if (std::abs(E - -1.110756115) > 1e-6){return 1;}
    else {return 0;}
}