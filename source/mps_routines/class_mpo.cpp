//
// Created by david on 2017-11-19.
//


#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_mpo.h>

using namespace std;
using namespace Textra;
using namespace Eigen;
using namespace std::complex_literals;

template<typename Scalar>
class_mpo<Scalar>::class_mpo() {
    local_dimension       = Model::local_dimension;
    h  = {Model::h(mps_sites,0), Model::h(mps_sites,1)};

    if constexpr (is_same_v<Scalar,std::complex<double>>){
        H_asMatrix      = Model::H(mps_sites);
        H_asTensor      = Textra::Matrix_to_Tensor<std::complex<double>,4>(Model::H(mps_sites), {2,2,2,2});
        H_MPO_asMatrix  = Model::MPO_asMatrix();
        H_MPO_asTensor2 = Textra::Matrix_to_Tensor2(Model::MPO_asMatrix());
    }else if constexpr (is_same_v<Scalar,double>){
        H_asMatrix      = Model::H(mps_sites).real();
        H_asTensor      = Textra::Matrix_to_Tensor<std::complex<double>,4>(Model::H(mps_sites), {2,2,2,2}).real();
        H_MPO_asMatrix  = Model::MPO_asMatrix().real();
        H_MPO_asTensor2 = Textra::Matrix_to_Tensor2(Model::MPO_asMatrix()).real();
    }else{
        std::cerr << "MPO not real or complex" << std::endl;
        exit(1);
    }

    Udt                   = compute_Udt(0.01, 1);
    M                     = compute_M();
    MM                    = compute_MM();
    F                     = compute_F(0.0001);
    G                     = compute_G(0.0001);


}



template <typename Scalar>
Tensor<Scalar,4> class_mpo<Scalar>::compute_M() {
    /*! Returns the MPO hamiltonian as a rank 4 MPO. Notation following Schollwöck (2010)
     *
     *          2
     *          |
     *      0---M---1
     *          |
     *          3
     */

    return Matrix_to_Tensor<Scalar,4> (H_MPO_asMatrix, {2,3,2,3}).shuffle(array4{1,3,0,2});

}

template <typename Scalar>
Tensor<Scalar,6> class_mpo<Scalar>::compute_MM() {
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

template <typename Scalar>
Tensor<std::complex<double>,4> class_mpo<Scalar>::compute_F(double a){
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


template <typename Scalar>
Tensor<std::complex<double>,4> class_mpo<Scalar>::compute_G(double a){
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


template <typename Scalar>
Tensor<std::complex<double>,4> class_mpo<Scalar>::compute_logG(double a){
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



template <typename Scalar>
Tensor<std::complex<double>,4> class_mpo<Scalar>::TimeEvolution_1st_order(const double delta_t) {
    return Matrix_to_Tensor<std::complex<double>,4>(Suzuki_Trotter_1st_order(-delta_t), array4{2,2,2,2});
}


template <typename Scalar>
Tensor<std::complex<double>,4> class_mpo<Scalar>::TimeEvolution_2nd_order(const double delta_t) {
    return Matrix_to_Tensor<std::complex<double>,4>(Suzuki_Trotter_2nd_order(-delta_t), array4{2,2,2,2});
}


template <typename Scalar>
Tensor<std::complex<double>,4> class_mpo<Scalar>::TimeEvolution_4th_order(const double delta_t) {
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

template <typename Scalar>
Tensor<Scalar,4> class_mpo<Scalar>::compute_Udt(double delta_t, int order) {
    /*! Returns a 2-site non-MPO time evolution operator.
*
*           0         1
*           |         |
*           [exp(H*dt)]
*           |         |
*           2         3
*/
    timestep = delta_t;
    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        switch (order) {
            case 1:
                return TimeEvolution_1st_order(delta_t);
            case 2:
                return TimeEvolution_2nd_order(delta_t);
            case 4:
                return TimeEvolution_4th_order(delta_t);
            default:
                return TimeEvolution_1st_order(delta_t);
        }
    } else if constexpr(std::is_same_v<Scalar, double>) {
        switch (order) {
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
}

template <typename Scalar>
void class_mpo<Scalar>::update_timestep(const double delta_t, const int susuki_trotter_order){
    timestep = delta_t;
    Udt = compute_Udt(delta_t, susuki_trotter_order);
}

//Explicit instantiations


template class class_mpo<double>;
template class class_mpo<std::complex<double>>;

