//
// Created by david on 2017-11-19.
//


#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <sim_parameters/nmspc_model.h>
#include <general/nmspc_tensor_extra.h>
#include "class_mpo.h"
#include <mps_routines/class_hamiltonian.h>
using namespace std;
using namespace Textra;
using namespace Eigen;
using namespace std::complex_literals;
using Scalar = class_mpo::Scalar;

class_mpo::class_mpo() {
//    local_dimension       = Model::local_dimension;
//    h                     = {Model::h(mps_sites,0), Model::h(mps_sites,1)};
//    H_asMatrix            = Model::H(mps_sites);
//    H_asTensor            = Textra::Matrix_to_Tensor(Model::H(mps_sites), 2,2,2,2);
//    H_asTensor_sq         = Textra::Matrix_to_Tensor(Model::Hsq(mps_sites), 2,2,2,2);


//    update_evolution_step_size(-0.01, 2);
//    HA                    = compute_H_MPO();
//    HB                    = compute_H_MPO();
//    G0                    = compute_G( 0.0i*1e-2, 4);
//    G1                    = compute_G( 1.0i*1e-2, 4);
}



//Tensor<Scalar,4> class_mpo::compute_H_MPO(double k)
///*! Returns the MPO hamiltonian as a rank 4 MPO. Notation following Schollwöck (2010)
// *
// *          2
// *          |
// *      0---H---1
// *          |
// *          3
// */
//{
//    return Matrix_to_Tensor(Model::H_MPO(k), 2,3,2,3).shuffle(array4{1,3,0,2});
//}
//
//Tensor<Scalar,4> class_mpo::compute_H_MPO_custom_field(double g, double e)
///*! Returns the MPO hamiltonian as a rank 4 MPO. Notation following Schollwöck (2010)
// *
// *          2
// *          |
// *      0---H---1
// *          |
// *          3
// */
//{
//    return Matrix_to_Tensor(Model::H_MPO_random_field(g,e), 2,3,2,3).shuffle(array4{1,3,0,2});
//}








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
