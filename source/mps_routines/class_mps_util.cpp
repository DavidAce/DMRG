//
// Created by david on 2018-07-06.
//

#include <mps_routines/class_mps_2site.h>
#include <general/class_eigsolver_arpack.h>
#include <general/nmspc_math.h>
#include "class_mps_util.h"

using namespace Textra;



void class_mps_util::compute_mps_components(const std::unique_ptr<class_mps_2site> &MPS){
//    int chiA2 = (int)(chiA()*chiA());
    int chiB2 = (int)(MPS->chiB()*MPS->chiB());
    int chiC2 = (int)(MPS->chiC()*MPS->chiC());


    Eigen::Tensor<Scalar,2> theta_evn_transfer_mat   = get_transfer_matrix_theta_evn(MPS).reshape(array2{chiB2,chiB2});
    Eigen::Tensor<Scalar,2> theta_odd_transfer_mat   = get_transfer_matrix_theta_odd(MPS).reshape(array2{chiC2,chiC2});
    class_eigsolver_arpack<Scalar, Form::GENERAL> solver;

    int ncvC = std::min(16, chiC2);
    int ncvB = std::min(16, chiB2);
    solver.eig(theta_evn_transfer_mat.data(), chiB2, 1, ncvB, Ritz::LM, Side::R, true, true);
    auto[eigvec_R_evn,eigval_R_evn]  = solver.get_eig_vecs_vals();
    solver.eig(theta_evn_transfer_mat.data(), chiB2, 1, ncvB, Ritz::LM, Side::L, true, true);
    auto[eigvec_L_evn,eigval_L_evn]  = solver.get_eig_vecs_vals();
    solver.eig(theta_odd_transfer_mat.data(), chiC2, 1, ncvC, Ritz::LM, Side::R, true, true);
    auto[eigvec_R_odd,eigval_R_odd]  = solver.get_eig_vecs_vals();
    solver.eig(theta_odd_transfer_mat.data(), chiC2, 1, ncvC, Ritz::LM, Side::L, true, true);
    auto[eigvec_L_odd,eigval_L_odd]  = solver.get_eig_vecs_vals();

    MatrixType<Scalar> eigvec_R_evn_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_R_evn.data(), eigvec_R_evn.size(),1);
    MatrixType<Scalar> eigvec_L_evn_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_L_evn.data(), eigvec_L_evn.size(),1);
    MatrixType<Scalar> eigvec_R_odd_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_R_odd.data(), eigvec_R_odd.size(),1);
    MatrixType<Scalar> eigvec_L_odd_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_L_odd.data(), eigvec_L_odd.size(),1);

    Scalar normalization_evn = sqrt((eigvec_L_evn_map.transpose() * eigvec_R_evn_map).sum());
    Scalar normalization_odd = sqrt((eigvec_L_odd_map.transpose() * eigvec_R_odd_map).sum());

    r_evn = Matrix_to_Tensor2(eigvec_R_evn_map).reshape(array2{MPS->chiB(),MPS->chiB()})/normalization_evn;
    l_evn = Matrix_to_Tensor2(eigvec_L_evn_map).reshape(array2{MPS->chiB(),MPS->chiB()})/normalization_evn;
    r_odd = Matrix_to_Tensor2(eigvec_R_odd_map).reshape(array2{MPS->chiC(),MPS->chiC()})/normalization_odd;
    l_odd = Matrix_to_Tensor2(eigvec_L_odd_map).reshape(array2{MPS->chiC(),MPS->chiC()})/normalization_odd;

    theta                = get_theta(MPS);
    theta_sw             = get_theta_swapped(MPS);
    theta_evn_normalized = get_theta_evn(MPS, sqrt(eigval_R_evn[0]));
    theta_odd_normalized = get_theta_odd(MPS, sqrt(eigval_R_odd[0]));

    LBGA                 = MPS->A();// / (T) sqrt(eigval_R_LBGA(0));
    LAGB                 = MPS->C().contract(MPS->MPS_B->get_G(), idx({1},{1})).shuffle(array3{1,0,2});// / (T) sqrt(eigval_R_LAGB(0));

    transfer_matrix_evn    = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    transfer_matrix_odd    = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<Scalar,4> transfer_matrix_LAGB_unnormalized = LAGB.contract(LAGB.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<Scalar,4> transfer_matrix_LBGA_unnormalized = LBGA.contract(LBGA.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});

    Eigen::Tensor<Scalar,0> l_evn_LBGA_r_odd = l_evn.contract(transfer_matrix_LBGA_unnormalized, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1}));
    Eigen::Tensor<Scalar,0> l_odd_LAGB_r_evn = l_odd.contract(transfer_matrix_LAGB_unnormalized, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1}));

    transfer_matrix_LAGB = transfer_matrix_LAGB_unnormalized /l_odd_LAGB_r_evn(0);
    transfer_matrix_LBGA = transfer_matrix_LBGA_unnormalized /l_evn_LBGA_r_odd(0);
    LAGB = LAGB / sqrt(l_odd_LAGB_r_evn(0));
    LBGA = LBGA / sqrt(l_evn_LBGA_r_odd(0));

//    std::cout << "Check:" << setprecision(10) <<  std::endl;
//    std::cout << " l_odd_LAGB_r_evn          = " << l_odd_LAGB_r_evn(0) << std::endl;
//    std::cout << " l_evn_LBGA_r_odd          = " << l_evn_LBGA_r_odd(0) << std::endl;
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




Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_theta(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LA ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return
            MPS->A().contract(MPS->C(), idx({2},{0}))
                    .contract(MPS->B(), idx({2},{1})) / norm;
}



Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_theta_swapped(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  MPS->C() //whatever L_A was in the previous step
            .contract(MPS->B(),            idx({1},{1}))
            .contract(MPS->MPS_A->get_G(), idx({2},{1}))
            .contract(MPS->C(), idx({3},{0}))
            .shuffle(array4{1,0,2,3})
            /norm;
}





Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_theta_evn(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return  MPS->A()
            .contract(MPS->C(),  idx({2},{0}))
            .contract(MPS->MPS_B->get_G(),  idx({2},{1}))
            //            .shuffle(array4{1,0,2,3})
            /norm;
}

Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_theta_odd(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]-- [ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  MPS->C()
            .contract(MPS->MPS_B->get_G(),         idx({1},{1}))
            .contract(MPS->A(),                    idx({2},{1}))
            .shuffle(array4{1,0,2,3})
            /norm;
}


Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_zero(const std::unique_ptr<class_mps_2site> &MPS) const {
    Eigen::Tensor<Scalar,1> I = MPS->LC;
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
}



Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_LBGA(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const {
    return MPS->A().contract( MPS->A().conjugate() , idx({0},{0}))
                   .shuffle(array4{0,3,1,2})
           /norm;
}


Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_GALC(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const {
    return MPS->C()
           .contract(MPS->MPS_A->get_G(),               idx({2},{0}))
           .contract(MPS->MPS_A->get_G().conjugate(),   idx({0},{0}))
           .contract(MPS->C(),                          idx({3},{0}) )
           .shuffle(array4{0,2,1,3})
           /norm;
}

Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_GBLB(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const {
    return MPS->B().contract(MPS->B().conjugate() ,   idx({0},{0}))
                   .shuffle(array4{0,2,1,3})
           /norm;
}


Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_LCGB(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const {
    return  MPS->C()
             .contract(MPS->MPS_B->get_G(),               idx({1},{1}))
             .contract(MPS->MPS_B->get_G().conjugate(),   idx({1},{0}))
             .contract(MPS->C(),                          idx({2},{1}) )
             .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_theta_evn(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const {
    return get_theta_evn(MPS).contract(get_theta_evn(MPS).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_theta_odd(const std::unique_ptr<class_mps_2site> &MPS, Scalar norm) const {
    return get_theta_odd(MPS).contract(get_theta_odd(MPS).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Eigen::Tensor<class_mps_util::Scalar,4> class_mps_util::get_transfer_matrix_AB(const std::unique_ptr<class_mps_2site> &MPS, int p) const {
    Eigen::Tensor<Scalar,4> temp = get_transfer_matrix_zero(MPS);
    Eigen::Tensor<Scalar,4> temp2;
    for (int i = 0; i < p-2; i++){
        if(Math::mod(i,2) == 0){
            temp2 = temp.contract(get_transfer_matrix_LBGA(MPS), idx({2,3},{0,1}));

        }else{
            temp2 = temp.contract(get_transfer_matrix_LCGB(MPS), idx({2,3},{0,1}));
        }
        temp = temp2;


    }
    return temp;
}
