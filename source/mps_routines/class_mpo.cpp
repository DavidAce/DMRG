//
// Created by david on 2017-11-19.
//

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include "general/n_tensor_extra.h"
#include "class_mpo.h"

using namespace std;
using namespace Textra;
using namespace Eigen;
using namespace std::complex_literals;


Tensor<class_mpo::Scalar,4> class_mpo::compute_M() {
    /*! Returns the MPO hamiltonian as a rank 4 MPO. Notation following Schollwöck (2010)
     *
     *          2
     *          |
     *      0---M---1
     *          |
     *          3
     */

    return Matrix_to_Tensor<double,4> (H_MPO_asMatrix, {2,3,2,3}).shuffle(array4{1,3,0,2});

}

Tensor<class_mpo::Scalar,6> class_mpo::compute_MM() {
    /*! Returns a 2-site Hamitlonian MPO of rank 6. Notation following Schollwöck (2010)
     *
     *           2   3
     *           |   |
     *       0---M---M---1
     *           |   |
     *           4   5
     */
//    cout << M.dimensions() << endl;

    return  M.contract(M, idx<1>({1},{0})).shuffle(array6{0,3,1,4,2,5});

}


Tensor<std::complex<double>,4> class_mpo::compute_F(double a){
    /*! Returns the moment generating function MPO hamiltonian as a rank 4 MPO.
    *   F := exp(aM), where a is a small parameter and M is an MPO.
    *
    *           0       1
    *           |       |
    *           [exp(aH)]
    *           |       |
    *           2       3
    *
    * Where H = H_even + H_odd, so exp(aH) can be Suzuki-Trotter decomposed just as in iTEBD
    */
    return Matrix_to_Tensor<std::complex<double>,4>(Suzuki_Trotter_2nd_order(a), {2,2,2,2});
};



Tensor<std::complex<double>,4> class_mpo::compute_G(double a){
    /*! Returns the characteristic function MPO hamiltonian as a rank 4 tensor.
    *   G := exp(iaM), where a is a small parameter and M is an MPO.
    *           0         1
    *           |         |
    *           [exp(iaH)]
    *           |         |
    *           2         3
    */
    return Matrix_to_Tensor<std::complex<double>,4>(Suzuki_Trotter_1st_order(1.0i * a), {2,2,2,2});
};

Tensor<std::complex<double>,4> class_mpo::compute_logG(double a){
    /*! Returns the characteristic function MPO hamiltonian as a rank 4 tensor.
    *   G := exp(iaM), where a is a small parameter and M is an MPO.
    *           0         1
    *           |         |
    *           [exp(iaH)]
    *           |         |
    *           2         3
    */
    return Matrix_to_Tensor<std::complex<double>,4>(((1i*a*h[0]/2.0).exp() * (1i*a*h[1]).exp() * (1i*a * h[0]/2.0).exp()).log(), {2,2,2,2});
};



Tensor<std::complex<double>,4> class_mpo::TimeEvolution_1st_order(const double delta_t) {
    return Matrix_to_Tensor<std::complex<double>,4>(Suzuki_Trotter_1st_order(-delta_t), array4{2,2,2,2});
}

Tensor<std::complex<double>,4> class_mpo::TimeEvolution_2nd_order(const double delta_t) {
    return Matrix_to_Tensor<std::complex<double>,4>(Suzuki_Trotter_2nd_order(-delta_t), array4{2,2,2,2});
}


Tensor<std::complex<double>,4> class_mpo::TimeEvolution_4th_order(const double delta_t) {
    double delta1 = delta_t /(4.0-pow(4.0,1.0/3.0));
    double delta2 = delta1;
    double delta3 = delta_t - 2.0*delta1 - 2.0*delta2;
    return Matrix_to_Tensor<std::complex<double>,4>( Suzuki_Trotter_2nd_order(-delta1)
                                                    *Suzuki_Trotter_2nd_order(-delta2)
                                                    *Suzuki_Trotter_2nd_order(-delta3)
                                                    *Suzuki_Trotter_2nd_order(-delta2)
                                                    *Suzuki_Trotter_2nd_order(-delta1),
                                                     array4{2,2,2,2});

}

Tensor<class_mpo::Scalar,4> class_mpo::compute_Udt(double delta_t, int order){
    /*! Returns a 2-site non-MPO time evolution operator.
*
*           0         1
*           |         |
*           [exp(H*dt)]
*           |         |
*           2         3
*/
    timestep = delta_t;

    switch (order){
        case 1:
            return TimeEvolution_1st_order(delta_t).real();
        case 2:
            return TimeEvolution_2nd_order(delta_t).real();
        case 4:
            return TimeEvolution_4th_order(delta_t).real();
        default:
            return TimeEvolution_1st_order(delta_t).real();
    }
};


void class_mpo::update_timestep(const double delta_t, const int order){
    timestep = delta_t;
    Udt = compute_Udt(delta_t, order);
}