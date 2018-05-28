//
// Created by david on 2017-11-13.
//


#include <mps_routines/class_mps.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/class_arpack_eigsolver.h>
#include <iomanip>

using namespace std;
using namespace Textra;

using Scalar = class_mps::Scalar;

// Function definitions

void class_mps::initialize(const long local_dimension_){
    local_dimension = local_dimension_;
    long d = local_dimension;
    // Default is a random product state
    GA.resize(array3{d,1,1});
    GA.setZero();
    LA.resize(array1{1});
    LA.setConstant(1.0);
    GB.resize(array3{d,1,1});
    GB.setZero();
    LB.resize(array1{1});
    LB.setConstant(1.0);
    //Initialize to as spinors
    auto r1 = rn::uniform_complex_1();
    auto r2 = rn::uniform_complex_1();
    GA(1, 0, 0) = r1.real();
    GA(0, 0, 0) = r1.imag();
    GB(1, 0, 0) = r2.real();
    GB(0, 0, 0) = r2.imag();
    swapped = false;
    LB_left = LB;
    theta = get_theta();
}

void class_mps::swap_AB(){
    swapped    = !swapped;

    //Swap Gamma
    tmp3    = GA;
    GA      = GB;
    GB      = tmp3;

    //Swap Lambda
    tmp1    = LA;
    LA      = LB;
    LB      = tmp1;
    LB_left  = LB;

}

void class_mps::compute_mps_components(){
//    class_arpackpp_wrapper eig;

    long sizeLA = LA.size();
    long sizeLB = LB.size();
//    Tensor<T,2> LBGA_mat            = get_transfer_matrix_LBGA().reshape(array2{sizeLB*sizeLB,sizeLA*sizeLA});
//    Tensor<T,2> LAGB_mat            = get_transfer_matrix_LAGB().reshape(array2{sizeLA*sizeLA,sizeLB*sizeLB});
    Tensor<Scalar,2> theta_evn_transfer_mat       = get_transfer_matrix_theta_evn()  .reshape(array2{sizeLB*sizeLB,sizeLB*sizeLB});
    Tensor<Scalar,2> theta_odd_transfer_mat       = get_transfer_matrix_theta_odd()  .reshape(array2{sizeLA*sizeLA,sizeLA*sizeLA});
//    auto eigval_R_LBGA                   = eig.solve_dominant<arpack::Form::COMPLEX_GENERAL, arpack::Ritz::LM, arpack::Side::R, false>(LBGA_mat.data(), sizeLB*sizeLB,sizeLA*sizeLA, 1);
//    auto eigval_R_LAGB                   = eig.solve_dominant<arpack::Form::COMPLEX_GENERAL, arpack::Ritz::LM, arpack::Side::R, false>(LAGB_mat.data(), sizeLA*sizeLA,sizeLB*sizeLB, 1);
//    auto [eigvec_R_evn,eigval_R_evn]     = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, true>(theta_evn_transfer_mat.data(),sizeLB*sizeLB,sizeLB*sizeLB, 1);
//    auto [eigvec_L_evn,eigval_L_evn]     = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::L, true>(theta_evn_transfer_mat.data(),sizeLB*sizeLB,sizeLB*sizeLB, 1);
//    auto [eigvec_R_odd,eigval_R_odd]     = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, true>(theta_odd_transfer_mat.data(),sizeLA*sizeLA,sizeLA*sizeLA, 1);
//    auto [eigvec_L_odd,eigval_L_odd]     = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::L, true>(theta_odd_transfer_mat.data(),sizeLA*sizeLA,sizeLA*sizeLA, 1);

    class_arpack_eigsolver<Scalar, Form::GENERAL> solver;

    int nA = (int)(sizeLA*sizeLA);
    int nB = (int)(sizeLB*sizeLB);
    int ncvA = std::min(10, nA);
    int ncvB = std::min(10, nB);
    auto[eigvec_R_evn,eigval_R_evn]  = solver.eig_get_vec_val(theta_evn_transfer_mat.data(), Ritz::LM, Side::R,
                                                              sizeLB * sizeLB, 1, ncvB, true);
    auto[eigvec_L_evn,eigval_L_evn]  = solver.eig_get_vec_val(theta_evn_transfer_mat.data(), Ritz::LM, Side::L,
                                                              sizeLB * sizeLB, 1, ncvB, true);
    auto[eigvec_R_odd,eigval_R_odd]  = solver.eig_get_vec_val(theta_odd_transfer_mat.data(), Ritz::LM, Side::R,
                                                              sizeLA * sizeLA, 1, ncvA, true);
    auto[eigvec_L_odd,eigval_L_odd]  = solver.eig_get_vec_val(theta_odd_transfer_mat.data(), Ritz::LM, Side::L,
                                                              sizeLA * sizeLA, 1, ncvA, true);

    MatrixType<Scalar> eigvec_R_evn_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_R_evn.data(), eigvec_R_evn.size(),1);
    MatrixType<Scalar> eigvec_L_evn_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_L_evn.data(), eigvec_L_evn.size(),1);
    MatrixType<Scalar> eigvec_R_odd_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_R_odd.data(), eigvec_R_odd.size(),1);
    MatrixType<Scalar> eigvec_L_odd_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_L_odd.data(), eigvec_L_odd.size(),1);

    Scalar normalization_evn = sqrt((eigvec_L_evn_map.transpose() * eigvec_R_evn_map).sum());
    Scalar normalization_odd = sqrt((eigvec_L_odd_map.transpose() * eigvec_R_odd_map).sum());


    r_evn = Matrix_to_Tensor2(eigvec_R_evn_map).reshape(array2{sizeLB,sizeLB})/normalization_evn;
    l_evn = Matrix_to_Tensor2(eigvec_L_evn_map).reshape(array2{sizeLB,sizeLB})/normalization_evn;
    r_odd = Matrix_to_Tensor2(eigvec_R_odd_map).reshape(array2{sizeLA,sizeLA})/normalization_odd;
    l_odd = Matrix_to_Tensor2(eigvec_L_odd_map).reshape(array2{sizeLA,sizeLA})/normalization_odd;

    theta                = get_theta();
    theta_sw             = get_theta_swapped();
    theta_evn_normalized = get_theta_evn(sqrt(eigval_R_evn[0]));
    theta_odd_normalized = get_theta_odd(sqrt(eigval_R_odd[0]));

    LBGA                 = asDiagonal(LB).contract(GA, idx({1},{1})).shuffle(array3{1,0,2});// / (T) sqrt(eigval_R_LBGA(0));
    LAGB                 = asDiagonal(LA).contract(GB, idx({1},{1})).shuffle(array3{1,0,2});// / (T) sqrt(eigval_R_LAGB(0));

    transfer_matrix_evn    = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    transfer_matrix_odd    = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    Tensor<Scalar,4> transfer_matrix_LAGB_unnormalized = LAGB.contract(LAGB.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});
    Tensor<Scalar,4> transfer_matrix_LBGA_unnormalized = LBGA.contract(LBGA.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});

    Tensor<Scalar,0> l_evn_LBGA_r_odd = l_evn.contract(transfer_matrix_LBGA_unnormalized, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1}));
    Tensor<Scalar,0> l_odd_LAGB_r_evn = l_odd.contract(transfer_matrix_LAGB_unnormalized, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1}));


    transfer_matrix_LAGB = transfer_matrix_LAGB_unnormalized /l_odd_LAGB_r_evn(0);
    transfer_matrix_LBGA = transfer_matrix_LBGA_unnormalized /l_evn_LBGA_r_odd(0);
    LAGB = LAGB / sqrt(l_odd_LAGB_r_evn(0));
    LBGA = LBGA / sqrt(l_evn_LBGA_r_odd(0));


//    std::cout << "Check:" << std::endl;
//    std::cout << " < l_evn | r_evn >         = " << l_evn.contract(r_evn, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < l_odd | r_odd >         = " << l_odd.contract(r_odd, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < l_evn | LBGA  | r_odd > = " << l_evn.contract(transfer_matrix_LBGA, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < l_odd | LAGB  | r_evn > = " << l_odd.contract(transfer_matrix_LAGB, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta     | theta >     = " << theta.contract(theta.conjugate(), idx({1,3,0,2},{1,3,0,2})) << std::endl;
//    std::cout << " < theta_evn_normalized | theta_evn_normalized > = " << theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).contract(l_evn, idx({0,2},{0,1})).contract(r_evn,idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta_evn_normalized | theta_evn_normalized > = " << transfer_matrix_evn.contract(l_evn, idx({0,1},{0,1})).contract(r_evn,idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta_odd_normalized | theta_odd_normalized > = " << theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).contract(l_odd, idx({0,2},{0,1})).contract(r_odd,idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta_odd_normalized | theta_odd_normalized > = " << transfer_matrix_odd.contract(l_odd, idx({0,1},{0,1})).contract(r_odd,idx({0,1},{0,1})) << std::endl;


}


//
//Tensor<T,4> class_mps::get_thetaL(T norm) const
///*!
// * Returns a left normalized two-site MPS
//     @verbatim
//        1--[ GA ]--[ C ]-- [ GB ] -- [ LB ]--3
//             |                 |
//             0                 2
//     @endverbatim
// */
//{
//    return GA
//            .contract(asDiagonal(C), idx({2},{0}))
//            .contract(GB,             idx({2},{1}))
//            .contract(asDiagonal(LB), idx({3},{0}))
//            /norm;
//};




Tensor<Scalar,4> class_mps::get_theta(Scalar norm) const
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LA ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  asDiagonal(LB_left) //whatever L_A was in the previous step
            .contract(GA,             idx({1},{1}))
            .contract(asDiagonal(LA), idx({2},{0}))
            .contract(GB,             idx({2},{1}))
            .contract(asDiagonal(LB), idx({3},{0}))
            .shuffle(array4{1,0,2,3})
            /norm;

}



Tensor<Scalar,4> class_mps::get_theta_swapped(Scalar norm) const
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  asDiagonal(LA) //whatever L_A was in the previous step
            .contract(GB,             idx({1},{1}))
            .contract(asDiagonal(LB), idx({2},{0}))
            .contract(GA,             idx({2},{1}))
            .contract(asDiagonal(LA), idx({3},{0}))
            .shuffle(array4{1,0,2,3})
            /norm;
}


Tensor<Scalar,3> class_mps::A() const{
    return asDiagonal(LB_left).contract(GA, idx({1},{1})).shuffle(array3{1,0,2});
};

Tensor<Scalar,3> class_mps::B() const{
    return GB.contract(asDiagonal(LB), idx({2},{0}));
};




Tensor<Scalar,4> class_mps::get_theta_evn(Scalar norm) const
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]-- [ LA ] -- [ GB ]--3
             |                 |
             0                 2
     @endverbatim
 */

{
    return  asDiagonal(LB) //whatever L_A was in the previous step
            .contract(GA,             idx({1},{1}))
            .contract(asDiagonal(LA), idx({2},{0}))
            .contract(GB,             idx({2},{1}))
            .shuffle(array4{1,0,2,3})
            /norm;
};

Tensor<Scalar,4> class_mps::get_theta_odd(Scalar norm) const
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]-- [ LB ] -- [ GA ]--3
             |                 |
             0                 2
     @endverbatim
 */
{
    return  asDiagonal(LA)
            .contract(GB,               idx({1},{1}))
            .contract(asDiagonal(LB),   idx({2},{0}))
            .contract(GA,               idx({2},{1}))
            .shuffle(array4{1,0,2,3})
            /norm;
}


//
//Tensor<T,4> class_mps::get_thetaR_swapped() const
///*!
// * Returns a two-site MPS with A and B swapped
//     @verbatim
//        1--[ GB ]--[ LB ]-- [ GA ] -- [ C ]--3
//             |                 |
//             0                 2
//     @endverbatim
// */
//{
//    return  GB.contract(asDiagonal(LB), idx({2},{0}))
//            .contract(GA,               idx({2},{1}))
//            .contract(asDiagonal(C),   idx({3},{0}));
////            .shuffle(array4{1,0,2,3});
//}

Tensor<Scalar,4> class_mps::get_transfer_matrix_zero() const {
    Textra::Tensor<Scalar,1> I = LA;
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
};



Tensor<Scalar,4> class_mps::get_transfer_matrix_LBGA(Scalar norm) const {
    return asDiagonal(LB)
            .contract(GA,               idx({1},{1}))
            .contract(GA.conjugate(),   idx({1},{0}))
            .contract(asDiagonal(LB),   idx({2},{1}) )
            .shuffle(array4{0,3,1,2})
            /norm;
}


Tensor<Scalar,4> class_mps::get_transfer_matrix_GALA(Scalar norm) const {
    return asDiagonal(LA)
            .contract(GA,               idx({2},{0}))
            .contract(GA.conjugate(),   idx({0},{0}))
            .contract(asDiagonal(LA),   idx({3},{0}) )
            .shuffle(array4{0,2,1,3})
           /norm;
}

Tensor<Scalar,4> class_mps::get_transfer_matrix_GBLB(Scalar norm) const {
    return asDiagonal(LB)
            .contract(GB,                  idx({2},{0}))
            .contract(GB.conjugate(),      idx({0},{0}))
            .contract(asDiagonal(LB),      idx({3},{0}) )
            .shuffle(array4{0,2,1,3})
           /norm;
}


Tensor<Scalar,4> class_mps::get_transfer_matrix_LAGB(Scalar norm) const {
    return  asDiagonal(LA)
            .contract(GB,               idx({1},{1}))
            .contract(GB.conjugate(),   idx({1},{0}))
            .contract(asDiagonal(LA),   idx({2},{1}) )
            .shuffle(array4{0,3,1,2})
            /norm;
}


Tensor<Scalar,4> class_mps::get_transfer_matrix_theta_evn(Scalar norm) const {
    return get_theta_evn().contract(get_theta_evn().conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Tensor<Scalar,4> class_mps::get_transfer_matrix_theta_odd(Scalar norm) const {
    return get_theta_odd().contract(get_theta_odd().conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Tensor<Scalar,4> class_mps::get_transfer_matrix_AB(int p) const {
    Tensor<Scalar,4> temp = get_transfer_matrix_zero();
    Tensor<Scalar,4> temp2;
    for (int i = 0; i < p-2; i++){
        if(Math::mod(i,2) == 0){
            temp2 = temp.contract(get_transfer_matrix_LBGA(), idx({2,3},{0,1}));

        }else{
            temp2 = temp.contract(get_transfer_matrix_LAGB(), idx({2,3},{0,1}));
        }
        temp = temp2;


    }
    return temp;
}


//
//Tensor<T,4> class_mps::get_regularization_fixpointA() const {
//    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
//    return asDiagonalSquared(LB).contract(asDiagonalSquared(C),pair);
//}
//
//Tensor<T,4> class_mps::get_regularization_fixpointB() const {
//    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
//    return asDiagonalSquared(C).contract(asDiagonalSquared(LB),pair);
//}
//
//
//Tensor<T,4> class_mps::get_transfer_matrix_LBGA_regularized() const {
//    return get_transfer_matrix_LBGA() - get_regularization_fixpointA();
//}
//
//Tensor<T,4> class_mps::get_transfer_matrix_GALA_regularized() const {
//    return get_transfer_matrix_GALA() - get_regularization_fixpointB();
//}
//
//Tensor<T,4> class_mps::get_transfer_matrix_GBLB_regularized() const {
//    return get_transfer_matrix_GBLB() - get_regularization_fixpointA();
//}
//
//Tensor<T,4> class_mps::get_transfer_matrix_LAGB_regularized() const {
//    return get_transfer_matrix_LAGB() - get_regularization_fixpointB();
//}
//
//Tensor<T,4> class_mps::get_transfer_matrix_AB_regularized() const {
////    return get_transfer_matrix_theta_evn() ;
//
//    return get_transfer_matrix_LBGA_regularized()
//            .contract(get_transfer_matrix_LAGB_regularized(), idx({2,3},{0,1}));
////    return get_transfer_matrix_theta_evn();
//}
//
//Tensor<T,4> class_mps::get_transfer_matrix_BA_regularized() const {
//    return get_transfer_matrix_LAGB_regularized().contract(get_transfer_matrix_LBGA_regularized(), idx({2,3},{0,1}));
//
////    return get_transfer_matrix_theta_odd();
////    return get_transfer_matrix_theta_odd() - get_regularization_fixpointB();
//}
//
//
//
//Tensor<T,4> class_mps::get_transfer_matrix_AB_regularized_term(int p) const {
//    Tensor<T,4> temp = get_transfer_matrix_zero() ;
//    Tensor<T,4> temp2;
//    for (int i = 0; i < p-2; i++){
//        if(Math::mod(i,2) == 0){
//            temp2 = temp.contract(get_transfer_matrix_LBGA_regularized(), idx({2,3},{0,1}));
//
//        }else{
//            temp2 = temp.contract(get_transfer_matrix_GBLB_regularized(), idx({2,3},{0,1}));
//        }
//
//        temp = temp2;
//    }
//    return temp;
//}
//
//
//Textra::Tensor<T,4> class_mps::get_transfer_matrix_regularized_inverseA()const{
//    long sizeL = LB.size();
//    long sizeR = C.size();
//    MatrixType<T> one_minus_regularized_transfer_matrix =
//            MatrixType<T>::Identity(sizeL*sizeL, sizeR*sizeR)
//          - Tensor_to_Matrix(get_transfer_matrix_LBGA_regularized(), sizeL*sizeL,sizeR*sizeR );
////    MatrixType<Scalar> pinv = Ereg.completeOrthogonalDecomposition().pseudoInverse();
//    //Return the pseudoinverse
//    return Matrix_to_Tensor( one_minus_regularized_transfer_matrix.completeOrthogonalDecomposition().pseudoInverse() ,sizeL,sizeL,sizeR,sizeR);
//};
//
//
//Textra::Tensor<T,4> class_mps::get_transfer_matrix_regularized_inverseB()const{
//    long sizeL = C.size();
//    long sizeR = LB.size();
//    MatrixType<T> one_minus_regularized_transfer_matrix =
//            MatrixType<T>::Identity(sizeL*sizeL, sizeR*sizeR)
//            - Tensor_to_Matrix(get_transfer_matrix_GBLB_regularized(), sizeL*sizeL,sizeR*sizeR );
////    MatrixType<Scalar> pinv = Ereg.completeOrthogonalDecomposition().pseudoInverse();
//    //Return the pseudoinverse
//    return Matrix_to_Tensor( one_minus_regularized_transfer_matrix.completeOrthogonalDecomposition().pseudoInverse() ,sizeL,sizeL,sizeR,sizeR);
//};
//
//
//
//
//Textra::Tensor<T,4> class_mps::get_transfer_matrix_regularized_inverseAB()const{
//    long sizeL = LB.size();
//    long sizeR = LB.size();
//    MatrixType<T> one_minus_regularized_transfer_matrix =
//            MatrixType<T>::Identity(sizeL*sizeL, sizeR*sizeR)
//          - Tensor_to_Matrix(get_transfer_matrix_AB_regularized(), sizeL*sizeL,sizeR*sizeR );
//    //Return the pseudoinverse
//    return Matrix_to_Tensor( one_minus_regularized_transfer_matrix.completeOrthogonalDecomposition().pseudoInverse() ,sizeL,sizeL,sizeR,sizeR);
//};
//
//Textra::Tensor<T,4> class_mps::get_transfer_matrix_regularized_inverseBA()const{
//    long sizeL = C.size();
//    long sizeR = C.size();
//    MatrixType<T> one_minus_regularized_transfer_matrix =
//            MatrixType<T>::Identity(sizeL*sizeL, sizeR*sizeR)
//          - Tensor_to_Matrix(get_transfer_matrix_BA_regularized(), sizeL*sizeL,sizeR*sizeR );
//    //Return the pseudoinverse
//    return Matrix_to_Tensor( one_minus_regularized_transfer_matrix.completeOrthogonalDecomposition().pseudoInverse() ,sizeL,sizeL,sizeR,sizeR);
//};
//
//

//Tensor<T,4> class_mps::get_transfer_matrix_L()const{
//    return A().contract(A().conjugate(), idx({0},{0}))
//              .shuffle(array4{0,2,1,3});
//};
//
//Tensor<T,4> class_mps::get_transfer_matrix_R()const{
//    return B().contract(B().conjugate(), idx({0},{0}))
//              .shuffle(array4{0,2,1,3});;
//};
//
//
//Tensor<T,4> class_mps::get_transfer_2_site_matrix_L()const{
//    return thetaL().contract(get_thetaL().conjugate(), idx({0,1},{0,1})).shuffle(array4{0,2,1,3});
//};
//
//
//Tensor<T,4> class_mps::get_transfer_2_site_matrix_R()const{
//    return theta_evn_normalized().contract(get_theta_evn().conjugate(), idx({0,1},{0,1})).shuffle(array4{0,2,1,3});
//};
