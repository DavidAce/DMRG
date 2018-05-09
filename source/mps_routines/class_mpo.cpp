//
// Created by david on 2017-11-19.
//


#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <sim_parameters/nmspc_model.h>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_hamiltonian.h>
using namespace std;
using namespace Textra;
using namespace Eigen;
using namespace std::complex_literals;
using Scalar = class_mpo::Scalar;

class_mpo::class_mpo() {
    local_dimension       = Model::local_dimension;
    h                     = {Model::h(mps_sites,0), Model::h(mps_sites,1)};
    H_asMatrix            = Model::H(mps_sites);
    H_asTensor            = Textra::Matrix_to_Tensor(Model::H(mps_sites), 2,2,2,2);
    H_asTensor_sq         = Textra::Matrix_to_Tensor(Model::Hsq(mps_sites), 2,2,2,2);
    H_MPO_asMatrix        = Model::H_MPO();
    H_MPO_asTensor2       = Textra::Matrix_to_Tensor2(Model::H_MPO());


    update_evolution_step_size(-0.01, 2);
    HA                    = compute_H_MPO();
    HB                    = compute_H_MPO();
    G0                    = compute_G( 0.0i*1e-2, 4);
    G1                    = compute_G( 1.0i*1e-2, 4);
}



Tensor<Scalar,4> class_mpo::compute_H_MPO(double k)
/*! Returns the MPO hamiltonian as a rank 4 MPO. Notation following Schollwöck (2010)
 *
 *          2
 *          |
 *      0---H---1
 *          |
 *          3
 */
{
    return Matrix_to_Tensor(Model::H_MPO(k), 2,3,2,3).shuffle(array4{1,3,0,2});
}

Tensor<Scalar,4> class_mpo::compute_H_MPO_custom_field(double g, double e)
/*! Returns the MPO hamiltonian as a rank 4 MPO. Notation following Schollwöck (2010)
 *
 *          2
 *          |
 *      0---H---1
 *          |
 *          3
 */
{
    return Matrix_to_Tensor(Model::H_MPO_random_field(g,e), 2,3,2,3).shuffle(array4{1,3,0,2});
}





std::vector<Textra::Tensor<Scalar,4>> class_mpo::compute_G(Scalar a, int susuki_trotter_order)
/*! Returns the moment generating function, or characteristic function (if a is imaginary) for the Hamiltonian as a rank 4 tensor.
*   G := exp(iaM) or exp(aM), where a is a small parameter and M is an MPO.
*   Note that G(-a) = G(a)* if  exp(iaM) !
*
@verbatim
            0         1
            |         |
            [ exp(aH) ]
            |         |
            2         3
@endverbatim
*/
{
    return get_2site_evolution_gates(a,susuki_trotter_order);
};



std::vector<Textra::MatrixType<Scalar>> class_mpo::Suzuki_Trotter_1st_order(Scalar t){
    return {(t*h[0]).exp(),
            (t*h[1]).exp() };
}

std::vector<Textra::MatrixType<Scalar>> class_mpo::Suzuki_Trotter_2nd_order(Scalar t){
    return {(t*h[0]/2.0).exp(),
            (t*h[1]).exp(),
            (t*h[0]/2.0).exp()};
}


std::vector<Textra::MatrixType<Scalar>> class_mpo::Suzuki_Trotter_4th_order(Scalar t)
/*!
 * Implementation based on
 * Janke, W., & Sauer, T. (1992).
 * Properties of higher-order Trotter formulas.
 * Physics Letters A, 165(3), 199–205.
 * https://doi.org/10.1016/0375-9601(92)90035-K
 * */
{
    double cbrt2 = pow(2.0,1.0/3.0);
    double beta1 = 1.0/(2.0 - cbrt2);
    double beta2 = - cbrt2 *beta1;
    double alph1 = 0.5*beta1;
    double alph2 = (1.0 - cbrt2)/2.0 * beta1;

    std::vector<Textra::MatrixType<Scalar>> temp;

    temp.emplace_back( (alph1 *  t*h[0]).exp() );
    temp.emplace_back( (beta1 *  t*h[1]).exp() );
    temp.emplace_back( (alph2 *  t*h[0]).exp() );
    temp.emplace_back( (beta2 *  t*h[1]).exp() );
    temp.emplace_back( (alph2 *  t*h[0]).exp() );
    temp.emplace_back( (beta1 *  t*h[1]).exp() );
    temp.emplace_back( (alph1 *  t*h[0]).exp() );
    return temp;
}


std::vector<Textra::Tensor<Scalar,4>> class_mpo::get_2site_evolution_gates(const Scalar t,int susuki_trotter_order)
/*! Returns a set of 2-site unitary gates, using Suzuki Trotter decomposition to order 1, 2 or 3.
 * These gates need to be applied to the MPS one at a time with a swap in between.
 */
{
    std::vector<Textra::MatrixType<Scalar>> matrix_vec;
    switch (susuki_trotter_order) {
        case 1:  matrix_vec = Suzuki_Trotter_1st_order(t);break;
        case 2:  matrix_vec = Suzuki_Trotter_2nd_order(t);break;
        case 4:  matrix_vec = Suzuki_Trotter_4th_order(t);break;
        default: matrix_vec = Suzuki_Trotter_2nd_order(t);break;
    }
    std::vector<Textra::Tensor<Scalar ,4>> tensor_vec;
    for(auto &m : matrix_vec){
        tensor_vec.emplace_back(Textra::Matrix_to_Tensor(m, 2,2,2,2));
    }
    return tensor_vec;
}


void class_mpo::update_evolution_step_size(const Scalar dt, const int susuki_trotter_order){
    /*! Returns a set of 2-site unitary gates for the time evolution operator. */
    step_size = std::abs(dt);
    U = get_2site_evolution_gates(dt, susuki_trotter_order);
}




//
//Tensor<T,6> class_mpo::compute_MM(T k)
///*! Returns a 2-site Hamitlonian MPO of rank 6. Notation following Schollwöck (2010)
// *
// *           2   3
// *           |   |
// *       0---HA---HA---1
// *           |   |
// *           4   5
// */
//{
//    auto HA = compute_H_MPO(k);
//    return  HA.contract(HA, idx({1},{0})).shuffle(array6{0,3,1,4,2,5});
//}

//
//Tensor<T,8> class_mpo::compute_MMMM(T k)
///*! Returns a 2-site Hamitlonian MPO of rank 6. Notation following Schollwöck (2010)
// *
// *           4   5
// *           |   |
// *       0---HA---HA---1
// *           |   |
// *       2---HA---HA---3
// *           |   |
// *           6   7
// */
//{
//    auto MM = compute_MM(k);
//    return  MM.contract(MM, idx({4,5},{2,3})).shuffle(array8{0,1,4,5,2,3,6,7});
////    return  M.contract(M, idx({2,3},{4,5})).shuffle(array8{0,1,4,5,2,3,6,7});
//}
