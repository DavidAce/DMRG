//
// Created by david on 2019-01-31.
//



//
// Created by david on 2018-07-06.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <state/class_finite_state.h>
#include <state/class_mps_2site.h>
#include <math/class_eigsolver.h>
#include <math/nmspc_math.h>
//#include "class_mps_util.h"

using namespace Textra;

namespace tools::common::views{
    Eigen::Tensor<std::complex<double>,4> theta                  = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> theta_evn_normalized   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> theta_odd_normalized   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> theta_sw               = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,3> LAGA                   = Eigen::Tensor<std::complex<double>,3> ();
    Eigen::Tensor<std::complex<double>,3> LCGB                   = Eigen::Tensor<std::complex<double>,3> ();
    Eigen::Tensor<std::complex<double>,2> l_evn                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,2> r_evn                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,2> l_odd                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,2> r_odd                  = Eigen::Tensor<std::complex<double>,2> ();
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LAGA   = Eigen::Tensor<std::complex<double>,4> ();
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LCGB   = Eigen::Tensor<std::complex<double>,4> ();
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

void tools::common::views::compute_mps_components(const class_infinite_state & state){
    if (components_computed)return;
    // On even thetas we have  chiA = chiB on the outer bond
    // On odd thetas we have chiC on the outer bond.
    int chiA2 = (int)(state.MPS->chiA()*state.MPS->chiA());
    int chiB2 = (int)(state.MPS->chiB()*state.MPS->chiB());
    int chiC2 = (int)(state.MPS->chiC()*state.MPS->chiC());

    theta = get_theta(state);
    Eigen::Tensor<std::complex<double>,2> theta_evn_transfer_mat   = get_transfer_matrix_theta_evn(state).reshape(array2{chiB2,chiB2});
    Eigen::Tensor<std::complex<double>,2> theta_odd_transfer_mat   = get_transfer_matrix_theta_odd(state).reshape(array2{chiC2,chiC2});

    using namespace eigutils::eigSetting;
    int ncvA = std::min(16, chiA2);
    int ncvC = std::min(16, chiC2);
    int ncvB = std::min(16, chiB2);
    [[maybe_unused]] auto [eigvec_R_evn, eigval_R_evn] = dominant_eig<Side::R>(theta_evn_transfer_mat, chiB2, ncvB);
    [[maybe_unused]] auto [eigvec_L_evn, eigval_L_evn] = dominant_eig<Side::L>(theta_evn_transfer_mat, chiB2, ncvB);
    [[maybe_unused]] auto [eigvec_R_odd, eigval_R_odd] = dominant_eig<Side::R>(theta_odd_transfer_mat, chiC2, ncvC);
    [[maybe_unused]] auto [eigvec_L_odd, eigval_L_odd] = dominant_eig<Side::L>(theta_odd_transfer_mat, chiC2, ncvC);

    std::complex<double> normalization_evn = sqrt((eigvec_L_evn.transpose() * eigvec_R_evn).sum());
    std::complex<double> normalization_odd = sqrt((eigvec_L_odd.transpose() * eigvec_R_odd).sum());

    r_evn = MatrixTensorMap(eigvec_R_evn).reshape(array2{state.MPS->chiB(),state.MPS->chiB()})/normalization_evn;
    l_evn = MatrixTensorMap(eigvec_L_evn).reshape(array2{state.MPS->chiB(),state.MPS->chiB()})/normalization_evn;
    r_odd = MatrixTensorMap(eigvec_R_odd).reshape(array2{state.MPS->chiC(),state.MPS->chiC()})/normalization_odd;
    l_odd = MatrixTensorMap(eigvec_L_odd).reshape(array2{state.MPS->chiC(),state.MPS->chiC()})/normalization_odd;

    theta                = get_theta(state);
    theta_sw             = get_theta_swapped(state);
    theta_evn_normalized = get_theta_evn(state, sqrt(eigval_R_evn));
    theta_odd_normalized = get_theta_odd(state, sqrt(eigval_R_odd));

    LAGA                 = state.MPS->A_bare(); // Modified! May contain some error? Check git history for comparison
    LCGB                 = state.MPS->LC().contract(state.MPS->GB(), Textra::idx({1},{1})).shuffle(array3{1,0,2}); // Modified! May contain some error? Check git history for comparison

    transfer_matrix_evn    = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    transfer_matrix_odd    = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LCGB_unnormalized = LCGB.contract(LCGB.conjugate(), idx({0}, {0})).shuffle(array4{0, 2, 1, 3});
    Eigen::Tensor<std::complex<double>,4> transfer_matrix_LAGA_unnormalized = LAGA.contract(LAGA.conjugate(), idx({0}, {0})).shuffle(array4{0, 2, 1, 3});

    Eigen::Tensor<std::complex<double>,0> l_evn_LBGA_r_odd = l_evn.contract(transfer_matrix_LAGA_unnormalized, idx({0, 1}, {0, 1})).contract(r_odd, idx({0, 1}, {0, 1}));
    Eigen::Tensor<std::complex<double>,0> l_odd_LCGB_r_evn = l_odd.contract(transfer_matrix_LCGB_unnormalized, idx({0, 1}, {0, 1})).contract(r_evn, idx({0, 1}, {0, 1}));

    transfer_matrix_LCGB = transfer_matrix_LCGB_unnormalized / l_odd_LCGB_r_evn(0);
    transfer_matrix_LAGA = transfer_matrix_LAGA_unnormalized / l_evn_LBGA_r_odd(0);
    LCGB = LCGB / sqrt(l_odd_LCGB_r_evn(0));
    LAGA = LAGA / sqrt(l_evn_LBGA_r_odd(0));
    components_computed = true;
    std::cout << "Check:" << std::setprecision(10) <<  std::endl;
    std::cout << " l_odd_LCGB_r_evn          = " << l_odd_LCGB_r_evn(0) << std::endl;
    std::cout << " l_evn_LBGA_r_odd          = " << l_evn_LBGA_r_odd(0) << std::endl;
    std::cout << " < l_evn | r_evn >         = " << l_evn.contract(r_evn, idx({0,1},{0,1})) << std::endl;
    std::cout << " < l_odd | r_odd >         = " << l_odd.contract(r_odd, idx({0,1},{0,1})) << std::endl;
    std::cout << " < l_evn | LAGA  | r_odd > = " << l_evn.contract(transfer_matrix_LAGA, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1})) << std::endl;
    std::cout << " < l_odd | LCGB  | r_evn > = " << l_odd.contract(transfer_matrix_LCGB, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1})) << std::endl;
    std::cout << " < theta     | theta >     = " << theta.contract(theta.conjugate(), idx({1,3,0,2},{1,3,0,2})) << std::endl;
    std::cout << " < theta_evn_normalized | theta_evn_normalized > = " << theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).contract(l_evn, idx({0,2},{0,1})).contract(r_evn,idx({0,1},{0,1})) << std::endl;
    std::cout << " < theta_evn_normalized | theta_evn_normalized > = " << transfer_matrix_evn.contract(l_evn, idx({0,1},{0,1})).contract(r_evn,idx({0,1},{0,1})) << std::endl;
    std::cout << " < theta_odd_normalized | theta_odd_normalized > = " << theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).contract(l_odd, idx({0,2},{0,1})).contract(r_odd,idx({0,1},{0,1})) << std::endl;
    std::cout << " < theta_odd_normalized | theta_odd_normalized > = " << transfer_matrix_odd.contract(l_odd, idx({0,1},{0,1})).contract(r_odd,idx({0,1},{0,1})) << std::endl;
}



Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta(const class_finite_state & state, std::complex<double> norm)
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return state.MPS_L.back().get_M().contract(state.MPS_R.front().get_M(), Textra::idx({2},{1})) /norm;
}



Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta(const class_infinite_state & state, std::complex<double> norm)
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
            state.MPS->A().contract(state.MPS->B(), Textra::idx({2},{1})) / norm;
}



Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta_swapped(const class_infinite_state & state, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]--[ LB ]-- [ GA ] -- [ LC ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{

    return state.MPS->LC() //whatever L_A was in the previous moves
                    .contract(state.MPS->B() , idx({1},{1}))
                    .contract(state.MPS->GA(), idx({2},{1}))
                    .contract(state.MPS->LC(), idx({3}, {0}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}





Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta_evn(const class_infinite_state & state, std::complex<double> norm)
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LA ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return  state.MPS->A()
             .contract(state.MPS->GB(),  idx({2},{1}))
            /norm;
}

Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta_odd(const class_infinite_state & state, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]-- [ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return state.MPS->LC()
                    .contract(state.MPS->B(),    idx({1},{1}))
                    .contract(state.MPS->GA(),   idx({2},{1}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_zero(const class_infinite_state & state) {
    Eigen::Tensor<std::complex<double>,1> I = state.MPS->MPS_A->get_LC();
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
}



Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_LBGA(const class_infinite_state & state, std::complex<double> norm)  {
    return state.MPS->A().contract( state.MPS->A().conjugate() , idx({0},{0}))
                   .shuffle(array4{0,3,1,2})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_GALC(const class_infinite_state & state, std::complex<double> norm)  {
    return state.MPS->LC()
                   .contract(state.MPS->GA(),             idx({2},{0}))
                   .contract(state.MPS->GA().conjugate(), idx({0},{0}))
                   .contract(state.MPS->LC(), idx({3}, {0}) )
                   .shuffle(array4{0,2,1,3})
           /norm;
}

Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_GBLB(const class_infinite_state & state, std::complex<double> norm)  {
    return state.MPS->B().contract(state.MPS->B().conjugate() ,   idx({0},{0}))
                   .shuffle(array4{0,2,1,3})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_LCGB(const class_infinite_state & state, std::complex<double> norm)  {
    return state.MPS->LC()
                    .contract(state.MPS->GB(),               idx({1},{1}))
                    .contract(state.MPS->GB().conjugate(),   idx({1},{0}))
                    .contract(state.MPS->LC(), idx({2}, {1}) )
                    .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_theta_evn(const class_infinite_state & state, std::complex<double> norm)  {
    using namespace tools::common::views;
    return get_theta_evn(state).contract(get_theta_evn(state).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_theta_odd(const class_infinite_state & state, std::complex<double> norm)  {
    return get_theta_odd(state).contract(get_theta_odd(state).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_AB(const class_infinite_state & state, int p) {
    Eigen::Tensor<std::complex<double>,4> temp = get_transfer_matrix_zero(state);
    Eigen::Tensor<std::complex<double>,4> temp2;
    for (int i = 0; i < p-2; i++){
        if(math::mod(i,2) == 0){
            temp2 = temp.contract(get_transfer_matrix_LBGA(state), idx({2,3},{0,1}));

        }else{
            temp2 = temp.contract(get_transfer_matrix_LCGB(state), idx({2,3},{0,1}));
        }
        temp = temp2;


    }
    return temp;
}



















Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LA ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return
            MPS.A().contract(MPS.B(), idx({2},{1})) / norm;
}



Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta_swapped(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return MPS.LC() //whatever L_A was in the previous moves
                    .contract(MPS.B(),        idx({1},{1}))
                    .contract(MPS.GA(),       idx({2},{1}))
                    .contract(MPS.LA(),       idx({3},{0}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}





Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta_evn(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LA ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return  MPS.A()
             .contract(MPS.GB(),  idx({2},{1}))
            /norm;
}

Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_theta_odd(const class_mps_2site  &MPS, std::complex<double> norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]-- [ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return MPS.LC()
                    .contract(MPS.B(),           idx({1},{1}))
                    .contract(MPS.GA(),          idx({2},{1}))
                    .shuffle(array4{1,0,2,3})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_zero(const class_mps_2site  &MPS) {
    Eigen::Tensor<std::complex<double>,1> I = MPS.MPS_A->get_LC();
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
}



Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_LBGA(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.A().contract(MPS.A().conjugate() , idx({0},{0}))
                   .shuffle(array4{0,3,1,2})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_GALC(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.LC()
                   .contract(MPS.GA(),               idx({2},{0}))
                   .contract(MPS.GA().conjugate(),   idx({0},{0}))
                   .contract(MPS.LC(), idx({3}, {0}) )
                   .shuffle(array4{0,2,1,3})
           /norm;
}

Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_GBLB(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.B().contract(MPS.B().conjugate() ,   idx({0},{0}))
                   .shuffle(array4{0,2,1,3})
           /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_LCGB(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return MPS.LC()
                    .contract(MPS.GB(),               idx({1},{1}))
                    .contract(MPS.GB().conjugate(),   idx({1},{0}))
                    .contract(MPS.LC(), idx({2}, {1}) )
                    .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_theta_evn(const class_mps_2site  &MPS, std::complex<double> norm)  {
    using namespace tools::common::views;
    return get_theta_evn(MPS).contract(get_theta_evn(MPS).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_theta_odd(const class_mps_2site  &MPS, std::complex<double> norm)  {
    return get_theta_odd(MPS).contract(get_theta_odd(MPS).conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Eigen::Tensor<std::complex<double>,4>
tools::common::views::get_transfer_matrix_AB(const class_mps_2site  &MPS, int p) {
    Eigen::Tensor<std::complex<double>,4> temp = get_transfer_matrix_zero(MPS);
    Eigen::Tensor<std::complex<double>,4> temp2;
    for (int i = 0; i < p-2; i++){
        if(math::mod(i,2) == 0){
            temp2 = temp.contract(get_transfer_matrix_LBGA(MPS), idx({2,3},{0,1}));

        }else{
            temp2 = temp.contract(get_transfer_matrix_LCGB(MPS), idx({2,3},{0,1}));
        }
        temp = temp2;


    }
    return temp;
}



