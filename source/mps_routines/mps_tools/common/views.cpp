//
// Created by david on 2019-01-31.
//



//
// Created by david on 2018-07-06.
//

#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_mps_2site.h>
#include <general/class_eigsolver.h>
#include <general/nmspc_math.h>
//#include "class_mps_util.h"

using namespace Textra;

namespace MPS_Tools::Common::Views{
    Eigen::Tensor<std::complex<double>,4> theta                  = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> theta_evn_normalized   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> theta_odd_normalized   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> theta_sw               = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,3> LBGA                   = Eigen::Tensor<std::complex<double>,3> ();
    Eigen::Tensor<std::complex<double>,3> LAGB                   = Eigen::Tensor<std::complex<double>,3> ();
    Eigen::Tensor<std::complex<double>,2> l_evn                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,2> r_evn                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,2> l_odd                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,2> r_odd                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LBGA   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LAGB   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_evn    = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_odd    = Eigen::Tensor<std::complex<double>,4> ();
    bool components_computed = false;
}

template<eigutils::eigSetting::Side side>
std::pair<Eigen::VectorXcd, std::complex<double>> dominant_eig(Eigen::Tensor<std::complex<double>,2> transfer_mat, int L, int ncv){
    using namespace eigutils::eigSetting;
    class_eigsolver solver;
    solver.eigs<Storage::DENSE>(transfer_mat.data(),L, 1, ncv,NAN,Form::NONSYMMETRIC,Ritz::LM,side, true,true);
    Eigen::VectorXcd eigvec = Eigen::Map<const Eigen::VectorXcd>(solver.solution.get_eigvecs<Type::CPLX,Form::NONSYMMETRIC, side>().data(), solver.solution.meta.rows,1);
    std::complex<double> eigval= solver.solution.get_eigvals<Form::NONSYMMETRIC>()[0];
    return std::make_pair(eigvec,eigval);
}

void MPS_Tools::Common::Views::compute_mps_components(const class_superblock & superblock){
//    int chiA2 = (int)(chiA()*chiA());
    if (components_computed)return;

    int chiB2 = (int)(superblock.MPS->chiB()*superblock.MPS->chiB());
    int chiC2 = (int)(superblock.MPS->chiC()*superblock.MPS->chiC());

    theta = get_theta(superblock);
    Eigen::Tensor<std::complex<double>,2> theta_evn_transfer_mat   = get_transfer_matrix_theta_evn(superblock).reshape(array2{chiB2,chiB2});
    Eigen::Tensor<std::complex<double>,2> theta_odd_transfer_mat   = get_transfer_matrix_theta_odd(superblock).reshape(array2{chiC2,chiC2});

    using namespace eigutils::eigSetting;
    int ncvC = std::min(16, chiC2);
    int ncvB = std::min(16, chiB2);
    [[maybe_unused]] auto [eigvec_R_evn, eigval_R_evn] = dominant_eig<Side::R>(theta_evn_transfer_mat, chiB2, ncvB);
    [[maybe_unused]] auto [eigvec_L_evn, eigval_L_evn] = dominant_eig<Side::L>(theta_evn_transfer_mat, chiB2, ncvB);
    [[maybe_unused]] auto [eigvec_R_odd, eigval_R_odd] = dominant_eig<Side::R>(theta_odd_transfer_mat, chiC2, ncvC);
    [[maybe_unused]] auto [eigvec_L_odd, eigval_L_odd] = dominant_eig<Side::L>(theta_odd_transfer_mat, chiC2, ncvC);

    std::complex<double> normalization_evn = sqrt((eigvec_L_evn.transpose() * eigvec_R_evn).sum());
    std::complex<double> normalization_odd = sqrt((eigvec_L_odd.transpose() * eigvec_R_odd).sum());

    r_evn = Matrix_to_Tensor2(eigvec_R_evn).reshape(array2{superblock.MPS->chiB(),superblock.MPS->chiB()})/normalization_evn;
    l_evn = Matrix_to_Tensor2(eigvec_L_evn).reshape(array2{superblock.MPS->chiB(),superblock.MPS->chiB()})/normalization_evn;
    r_odd = Matrix_to_Tensor2(eigvec_R_odd).reshape(array2{superblock.MPS->chiC(),superblock.MPS->chiC()})/normalization_odd;
    l_odd = Matrix_to_Tensor2(eigvec_L_odd).reshape(array2{superblock.MPS->chiC(),superblock.MPS->chiC()})/normalization_odd;

    theta                = get_theta(superblock);
    theta_sw             = get_theta_swapped(superblock);
    theta_evn_normalized = get_theta_evn(superblock, sqrt(eigval_R_evn));
    theta_odd_normalized = get_theta_odd(superblock, sqrt(eigval_R_odd));

    LBGA                 = superblock.MPS->A();// / (Scalar_) sqrt(eigval_R_LBGA(0));
    LAGB                 = superblock.MPS->C().contract(superblock.MPS->MPS_B->get_G(), idx({1},{1})).shuffle(array3{1,0,2});// / (Scalar_) sqrt(eigval_R_LAGB(0));

    transfer_matrix_evn    = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    transfer_matrix_odd    = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LAGB_unnormalized = LAGB.contract(LAGB.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LBGA_unnormalized = LBGA.contract(LBGA.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});

    Eigen::Tensor<std::complex<double>,0> l_evn_LBGA_r_odd = l_evn.contract(transfer_matrix_LBGA_unnormalized, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1}));
    Eigen::Tensor<std::complex<double>,0> l_odd_LAGB_r_evn = l_odd.contract(transfer_matrix_LAGB_unnormalized, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1}));

    transfer_matrix_LAGB = transfer_matrix_LAGB_unnormalized /l_odd_LAGB_r_evn(0);
    transfer_matrix_LBGA = transfer_matrix_LBGA_unnormalized /l_evn_LBGA_r_odd(0);
    LAGB = LAGB / sqrt(l_odd_LAGB_r_evn(0));
    LBGA = LBGA / sqrt(l_evn_LBGA_r_odd(0));
    components_computed = true;
//    std::cout << "Check:" << std::setprecision(10) <<  std::endl;
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



Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta(const class_finite_chain_state & state, std::complex<double> norm)
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return
            state.MPS_L.back().get_A().contract(Textra::asDiagonal(state.MPS_C), idx({2},{0}))
                    .contract(state.MPS_R.front().get_B(), idx({2},{1})) / norm;
}



Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta(const class_superblock & superblock, std::complex<double> norm)
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
            superblock.MPS->A().contract(superblock.MPS->C(), idx({2},{0}))
                    .contract(superblock.MPS->B(), idx({2},{1})) / norm;
}



Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta_swapped(const class_superblock & superblock, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  superblock.MPS->C() //whatever L_A was in the previous step
                    .contract(superblock.MPS->B(),            idx({1},{1}))
                    .contract(superblock.MPS->MPS_A->get_G(), idx({2},{1}))
                    .contract(superblock.MPS->C(), idx({3},{0}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}





Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta_evn(const class_superblock & superblock, std::complex<double> norm)
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return  superblock.MPS->A()
                    .contract(superblock.MPS->C(),  idx({2},{0}))
                    .contract(superblock.MPS->MPS_B->get_G(),  idx({2},{1}))
            //            .shuffle(array4{1,0,2,3})
            /norm;
}

Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta_odd(const class_superblock & superblock, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]-- [ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  superblock.MPS->C()
                    .contract(superblock.MPS->MPS_B->get_G(),         idx({1},{1}))
                    .contract(superblock.MPS->A(),                    idx({2},{1}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_zero(const class_superblock & superblock) {
    Eigen::Tensor<std::complex<double>,1> I = superblock.MPS->LC;
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
}



Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_LBGA(const class_superblock & superblock, std::complex<double> norm)  {
    return superblock.MPS->A().contract( superblock.MPS->A().conjugate() , idx({0},{0}))
                   .shuffle(array4{0,3,1,2})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_GALC(const class_superblock & superblock, std::complex<double> norm)  {
    return superblock.MPS->C()
                   .contract(superblock.MPS->MPS_A->get_G(),               idx({2},{0}))
                   .contract(superblock.MPS->MPS_A->get_G().conjugate(),   idx({0},{0}))
                   .contract(superblock.MPS->C(),                          idx({3},{0}) )
                   .shuffle(array4{0,2,1,3})
           /norm;
}

Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_GBLB(const class_superblock & superblock, std::complex<double> norm)  {
    return superblock.MPS->B().contract(superblock.MPS->B().conjugate() ,   idx({0},{0}))
                   .shuffle(array4{0,2,1,3})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_LCGB(const class_superblock & superblock, std::complex<double> norm)  {
    return  superblock.MPS->C()
                    .contract(superblock.MPS->MPS_B->get_G(),               idx({1},{1}))
                    .contract(superblock.MPS->MPS_B->get_G().conjugate(),   idx({1},{0}))
                    .contract(superblock.MPS->C(),                          idx({2},{1}) )
                    .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_theta_evn(const class_superblock & superblock, std::complex<double> norm)  {
    using namespace MPS_Tools::Common::Views;
    return get_theta_evn(superblock).contract(get_theta_evn(superblock).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_theta_odd(const class_superblock & superblock, std::complex<double> norm)  {
    return get_theta_odd(superblock).contract(get_theta_odd(superblock).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_AB(const class_superblock & superblock, int p) {
    Eigen::Tensor<std::complex<double>,4> temp = get_transfer_matrix_zero(superblock);
    Eigen::Tensor<std::complex<double>,4> temp2;
    for (int i = 0; i < p-2; i++){
        if(Math::mod(i,2) == 0){
            temp2 = temp.contract(get_transfer_matrix_LBGA(superblock), idx({2,3},{0,1}));

        }else{
            temp2 = temp.contract(get_transfer_matrix_LCGB(superblock), idx({2,3},{0,1}));
        }
        temp = temp2;


    }
    return temp;
}



















Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta(const class_mps_2site  &MPS, std::complex<double> norm)
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
            MPS.A().contract(MPS.C(), idx({2},{0}))
                    .contract(MPS.B(), idx({2},{1})) / norm;
}



Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta_swapped(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  MPS.C() //whatever L_A was in the previous step
                    .contract(MPS.B(),            idx({1},{1}))
                    .contract(MPS.MPS_A->get_G(), idx({2},{1}))
                    .contract(MPS.C(), idx({3},{0}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}





Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta_evn(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return  MPS.A()
                    .contract(MPS.C(),  idx({2},{0}))
                    .contract(MPS.MPS_B->get_G(),  idx({2},{1}))
            //            .shuffle(array4{1,0,2,3})
            /norm;
}

Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_theta_odd(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]-- [ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  MPS.C()
                    .contract(MPS.MPS_B->get_G(),         idx({1},{1}))
                    .contract(MPS.A(),                    idx({2},{1}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_zero(const class_mps_2site  &MPS) {
    Eigen::Tensor<std::complex<double>,1> I = MPS.LC;
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
}



Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_LBGA(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.A().contract(MPS.A().conjugate() , idx({0},{0}))
                   .shuffle(array4{0,3,1,2})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_GALC(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.C()
                   .contract(MPS.MPS_A->get_G(),               idx({2},{0}))
                   .contract(MPS.MPS_A->get_G().conjugate(),   idx({0},{0}))
                   .contract(MPS.C(),                          idx({3},{0}) )
                   .shuffle(array4{0,2,1,3})
           /norm;
}

Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_GBLB(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.B().contract(MPS.B().conjugate() ,   idx({0},{0}))
                   .shuffle(array4{0,2,1,3})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_LCGB(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return  MPS.C()
                    .contract(MPS.MPS_B->get_G(),               idx({1},{1}))
                    .contract(MPS.MPS_B->get_G().conjugate(),   idx({1},{0}))
                    .contract(MPS.C(),                          idx({2},{1}) )
                    .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_theta_evn(const class_mps_2site  &MPS, std::complex<double> norm)  {
    using namespace MPS_Tools::Common::Views;
    return get_theta_evn(MPS).contract(get_theta_evn(MPS).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_theta_odd(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return get_theta_odd(MPS).contract(get_theta_odd(MPS).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Eigen::Tensor<std::complex<double>,4>
MPS_Tools::Common::Views::get_transfer_matrix_AB(const class_mps_2site  &MPS, int p) {
    Eigen::Tensor<std::complex<double>,4> temp = get_transfer_matrix_zero(MPS);
    Eigen::Tensor<std::complex<double>,4> temp2;
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



