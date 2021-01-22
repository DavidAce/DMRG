//
// Created by david on 2019-01-31.
//

//
// Created by david on 2018-07-06.
//

#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/num.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/views.h>

// using namespace Textra;
using Scalar = tools::common::views::Scalar;
namespace tools::common::views {
    Eigen::Tensor<Scalar, 4> theta                = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 4> theta_evn_normalized = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 4> theta_odd_normalized = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 4> theta_sw             = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 3> LAGA                 = Eigen::Tensor<Scalar, 3>();
    Eigen::Tensor<Scalar, 3> LCGB                 = Eigen::Tensor<Scalar, 3>();
    Eigen::Tensor<Scalar, 2> l_evn                = Eigen::Tensor<Scalar, 2>();
    Eigen::Tensor<Scalar, 2> r_evn                = Eigen::Tensor<Scalar, 2>();
    Eigen::Tensor<Scalar, 2> l_odd                = Eigen::Tensor<Scalar, 2>();
    Eigen::Tensor<Scalar, 2> r_odd                = Eigen::Tensor<Scalar, 2>();
    Eigen::Tensor<Scalar, 4> transfer_matrix_LAGA = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 4> transfer_matrix_LCGB = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 4> transfer_matrix_evn  = Eigen::Tensor<Scalar, 4>();
    Eigen::Tensor<Scalar, 4> transfer_matrix_odd  = Eigen::Tensor<Scalar, 4>();
    bool                     components_computed  = false;
}

template<eig::Side side>
std::pair<Eigen::VectorXcd, Scalar> dominant_eig(const Eigen::Tensor<Scalar, 2> &transfer_mat, int L, int ncv) {
    eig::solver solver;
    solver.eigs(transfer_mat.data(), L, 1, ncv, eig::Ritz::LM, eig::Form::NSYM, side, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF, eig::Dephase::ON);
    auto eigvec = eig::view::get_eigvec<Scalar>(solver.result, 0, side);
    auto eigval = eig::view::get_eigval<Scalar>(solver.result, 0);
    return std::make_pair(eigvec, eigval);
}

void tools::common::views::compute_mps_components(const class_state_infinite &state) {
    if(components_computed) return;
    // On even thetas we have  chiA = chiB on the outer bond
    // On odd thetas we have chiC on the outer bond.
    //    int chiA2 = (int)(state.MPS->chiA()*state.MPS->chiA());
    int chiB2 = static_cast<int>(state.chiB() * state.chiB());
    int chiC2 = static_cast<int>(state.chiC() * state.chiC());

    theta                                           = get_theta(state);
    Eigen::Tensor<Scalar, 2> theta_evn_transfer_mat = get_transfer_matrix_theta_evn(state).reshape(Textra::array2{chiB2, chiB2});
    Eigen::Tensor<Scalar, 2> theta_odd_transfer_mat = get_transfer_matrix_theta_odd(state).reshape(Textra::array2{chiC2, chiC2});

    //    int ncvA = std::min(16, chiA2);
    int ncvC                                           = std::min(16, chiC2);
    int ncvB                                           = std::min(16, chiB2);
    [[maybe_unused]] auto [eigvec_R_evn, eigval_R_evn] = dominant_eig<eig::Side::R>(theta_evn_transfer_mat, chiB2, ncvB);
    [[maybe_unused]] auto [eigvec_L_evn, eigval_L_evn] = dominant_eig<eig::Side::L>(theta_evn_transfer_mat, chiB2, ncvB);
    [[maybe_unused]] auto [eigvec_R_odd, eigval_R_odd] = dominant_eig<eig::Side::R>(theta_odd_transfer_mat, chiC2, ncvC);
    [[maybe_unused]] auto [eigvec_L_odd, eigval_L_odd] = dominant_eig<eig::Side::L>(theta_odd_transfer_mat, chiC2, ncvC);

    Scalar normalization_evn = sqrt((eigvec_L_evn.transpose() * eigvec_R_evn).sum());
    Scalar normalization_odd = sqrt((eigvec_L_odd.transpose() * eigvec_R_odd).sum());

    r_evn = Textra::MatrixTensorMap(eigvec_R_evn).reshape(Textra::array2{state.chiB(), state.chiB()}) / normalization_evn;
    l_evn = Textra::MatrixTensorMap(eigvec_L_evn).reshape(Textra::array2{state.chiB(), state.chiB()}) / normalization_evn;
    r_odd = Textra::MatrixTensorMap(eigvec_R_odd).reshape(Textra::array2{state.chiC(), state.chiC()}) / normalization_odd;
    l_odd = Textra::MatrixTensorMap(eigvec_L_odd).reshape(Textra::array2{state.chiC(), state.chiC()}) / normalization_odd;

    theta                = get_theta(state);
    theta_sw             = get_theta_swapped(state);
    theta_evn_normalized = get_theta_evn(state, sqrt(eigval_R_evn));
    theta_odd_normalized = get_theta_odd(state, sqrt(eigval_R_odd));

    LAGA = state.A_bare(); // Modified! May contain some error? Check git history for comparison
    LCGB = state.LC_diag()
               .contract(state.GB(), Textra::idx({1}, {1}))
               .shuffle(Textra::array3{1, 0, 2}); // Modified! May contain some error? Check git history for comparison

    transfer_matrix_evn = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), Textra::idx({0, 2}, {0, 2})).shuffle(Textra::array4{0, 2, 1, 3});
    transfer_matrix_odd = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), Textra::idx({0, 2}, {0, 2})).shuffle(Textra::array4{0, 2, 1, 3});
    Eigen::Tensor<Scalar, 4> transfer_matrix_LCGB_unnormalized = LCGB.contract(LCGB.conjugate(), Textra::idx({0}, {0})).shuffle(Textra::array4{0, 2, 1, 3});
    Eigen::Tensor<Scalar, 4> transfer_matrix_LAGA_unnormalized = LAGA.contract(LAGA.conjugate(), Textra::idx({0}, {0})).shuffle(Textra::array4{0, 2, 1, 3});

    Eigen::Tensor<Scalar, 0> l_evn_LBGA_r_odd =
        l_evn.contract(transfer_matrix_LAGA_unnormalized, Textra::idx({0, 1}, {0, 1})).contract(r_odd, Textra::idx({0, 1}, {0, 1}));
    Eigen::Tensor<Scalar, 0> l_odd_LCGB_r_evn =
        l_odd.contract(transfer_matrix_LCGB_unnormalized, Textra::idx({0, 1}, {0, 1})).contract(r_evn, Textra::idx({0, 1}, {0, 1}));

    transfer_matrix_LCGB = transfer_matrix_LCGB_unnormalized / l_odd_LCGB_r_evn(0);
    transfer_matrix_LAGA = transfer_matrix_LAGA_unnormalized / l_evn_LBGA_r_odd(0);
    LCGB                 = LCGB / sqrt(l_odd_LCGB_r_evn(0));
    LAGA                 = LAGA / sqrt(l_evn_LBGA_r_odd(0));
    components_computed  = true;
    if(tools::log->level() == spdlog::level::trace) {
        Eigen::Tensor<Scalar, 0> le_re     = l_evn.contract(r_evn, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> lo_ro     = l_odd.contract(r_odd, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> le_laga   = l_evn.contract(transfer_matrix_LAGA, Textra::idx({0, 1}, {0, 1})).contract(r_odd, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> lo_lcgb   = l_odd.contract(transfer_matrix_LCGB, Textra::idx({0, 1}, {0, 1})).contract(r_evn, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> thth      = theta.contract(theta.conjugate(), Textra::idx({1, 3, 0, 2}, {1, 3, 0, 2}));
        Eigen::Tensor<Scalar, 0> then_then = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), Textra::idx({0, 2}, {0, 2}))
                                                 .contract(l_evn, Textra::idx({0, 2}, {0, 1}))
                                                 .contract(r_evn, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> thon_thon = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), Textra::idx({0, 2}, {0, 2}))
                                                 .contract(l_odd, Textra::idx({0, 2}, {0, 1}))
                                                 .contract(r_odd, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> tren_le = transfer_matrix_evn.contract(l_evn, Textra::idx({0, 1}, {0, 1})).contract(r_evn, Textra::idx({0, 1}, {0, 1}));
        Eigen::Tensor<Scalar, 0> trod_lo = transfer_matrix_odd.contract(l_odd, Textra::idx({0, 1}, {0, 1})).contract(r_odd, Textra::idx({0, 1}, {0, 1}));
        std::string              str     = "Check\n";
        str +=
            fmt::format(" l_odd_LCGB_r_evn                                = {:.15f}{:+.15f}\n", std::real(l_odd_LCGB_r_evn(0)), std::imag(l_odd_LCGB_r_evn(0)));
        str +=
            fmt::format(" l_evn_LBGA_r_odd                                = {:.15f}{:+.15f}\n", std::real(l_evn_LBGA_r_odd(0)), std::imag(l_evn_LBGA_r_odd(0)));
        str += fmt::format(" < l_evn | r_evn >                               = {:.15f}{:+.15f}\n", std::real(le_re(0)),
                           std::imag(le_re(0))); //  l_evn.contract(r_evn, Textra::idx({0,1},{0,1})));
        str += fmt::format(" < l_odd | r_odd >                               = {:.15f}{:+.15f}\n", std::real(lo_ro(0)),
                           std::imag(lo_ro(0))); //  l_odd.contract(r_odd, Textra::idx({0,1},{0,1})));
        str +=
            fmt::format(" < l_evn | LAGA  | r_odd >                       = {:.15f}{:+.15f}\n", std::real(le_laga(0)),
                        std::imag(le_laga(0))); //  l_evn.contract(transfer_matrix_LAGA, Textra::idx({0,1},{0,1})).contract(r_odd, Textra::idx({0,1},{0,1})));
        str +=
            fmt::format(" < l_odd | LCGB  | r_evn >                       = {:.15f}{:+.15f}\n", std::real(lo_lcgb(0)),
                        std::imag(lo_lcgb(0))); //  l_odd.contract(transfer_matrix_LCGB, Textra::idx({0,1},{0,1})).contract(r_evn, Textra::idx({0,1},{0,1})));
        str += fmt::format(" < theta     | theta >                           = {:.15f}{:+.15f}\n", std::real(thth(0)),
                           std::imag(thth(0))); //  theta.contract(theta.conjugate(), Textra::idx({1,3,0,2},{1,3,0,2})));
        str +=
            fmt::format(" < theta_evn_normalized | theta_evn_normalized > = {:.15f}{:+.15f}\n", std::real(then_then(0)),
                        std::imag(then_then(0))); //  theta_evn_normalized.contract(theta_evn_normalized.conjugate(), Textra::idx({0,2},{0,2})).contract(l_evn,
                                                  //  Textra::idx({0,2},{0,1})).contract(r_evn,Textra::idx({0,1},{0,1})));
        str +=
            fmt::format(" < theta_odd_normalized | theta_odd_normalized > = {:.15f}{:+.15f}\n", std::real(thon_thon(0)),
                        std::imag(thon_thon(0))); //  theta_odd_normalized.contract(theta_odd_normalized.conjugate(), Textra::idx({0,2},{0,2})).contract(l_odd,
                                                  //  Textra::idx({0,2},{0,1})).contract(r_odd,Textra::idx({0,1},{0,1})));
        str += fmt::format(" < theta_evn_normalized | theta_evn_normalized > = {:.15f}{:+.15f}\n", std::real(tren_le(0)),
                           std::imag(tren_le(0))); //  transfer_matrix_evn.contract(l_evn, Textra::idx({0,1},{0,1})).contract(r_evn,Textra::idx({0,1},{0,1})));
        str += fmt::format(" < theta_odd_normalized | theta_odd_normalized > = {:.15f}{:+.15f}\n", std::real(trod_lo(0)),
                           std::imag(trod_lo(0))); //  transfer_matrix_odd.contract(l_odd, Textra::idx({0,1},{0,1})).contract(r_odd,Textra::idx({0,1},{0,1})));
        tools::log->trace(str);
    }
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_theta(const class_state_finite &state, Scalar norm)
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2

    where we define LB*GA*LC = M on the "A" site and GB*LB = M on the right site
     @endverbatim
 */
{
    auto        pos  = state.get_position();
    const auto &mpsA = state.get_mps_site(pos);
    const auto &mpsB = state.get_mps_site(pos + 1);
    return mpsA.get_M().contract(mpsB.get_M(), Textra::idx({2}, {1})) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_theta(const class_state_infinite &state, Scalar norm)
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LA ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return state.A().contract(state.B(), Textra::idx({2}, {1})) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_theta_swapped(const class_state_infinite &state, Scalar norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]--[ LB ]-- [ GA ] -- [ LC ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return state
               .LC_diag() // whatever L_A was in the previous moves
               .contract(state.B(), Textra::idx({1}, {1}))
               .contract(state.GA(), Textra::idx({2}, {1}))
               .contract(state.LC_diag(), Textra::idx({3}, {0}))
               .shuffle(Textra::array4{1, 0, 2, 3}) /
           norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_theta_evn(const class_state_infinite &state, Scalar norm)
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LA ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return state.A().contract(state.GB(), Textra::idx({2}, {1})) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_theta_odd(const class_state_infinite &state, Scalar norm)
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LC ]--[ GB ]--[ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return state.LC_diag().contract(state.B(), Textra::idx({1}, {1})).contract(state.GA(), Textra::idx({2}, {1})).shuffle(Textra::array4{1, 0, 2, 3}) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_zero(const class_state_infinite &state) {
    long                     dim = state.chiC();
    Eigen::Tensor<Scalar, 2> I   = Textra::asDiagonal(Eigen::Tensor<Scalar, 1>(dim));
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>, 0> pair = {};
    return I.contract(I, pair).shuffle(Textra::array4{0, 2, 1, 3});
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_LBGA(const class_state_infinite &state, Scalar norm) {
    return state.A().contract(state.A().conjugate(), Textra::idx({0}, {0})).shuffle(Textra::array4{0, 3, 1, 2}) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_GALC(const class_state_infinite &state, Scalar norm) {
    return state.LC_diag()
               .contract(state.GA(), Textra::idx({2}, {0}))
               .contract(state.GA().conjugate(), Textra::idx({0}, {0}))
               .contract(state.LC_diag(), Textra::idx({3}, {0}))
               .shuffle(Textra::array4{0, 2, 1, 3}) /
           norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_GBLB(const class_state_infinite &state, Scalar norm) {
    return state.B().contract(state.B().conjugate(), Textra::idx({0}, {0})).shuffle(Textra::array4{0, 2, 1, 3}) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_LCGB(const class_state_infinite &state, Scalar norm) {
    return state.LC_diag()
               .contract(state.GB(), Textra::idx({1}, {1}))
               .contract(state.GB().conjugate(), Textra::idx({1}, {0}))
               .contract(state.LC_diag(), Textra::idx({2}, {1}))
               .shuffle(Textra::array4{0, 3, 1, 2}) /
           norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_theta_evn(const class_state_infinite &state, Scalar norm) {
    using namespace tools::common::views;
    return get_theta_evn(state).contract(get_theta_evn(state).conjugate(), Textra::idx({0, 2}, {0, 2})).shuffle(Textra::array4{0, 2, 1, 3}) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_theta_odd(const class_state_infinite &state, Scalar norm) {
    return get_theta_odd(state).contract(get_theta_odd(state).conjugate(), Textra::idx({0, 2}, {0, 2})).shuffle(Textra::array4{0, 2, 1, 3}) / norm;
}

Eigen::Tensor<Scalar, 4> tools::common::views::get_transfer_matrix_AB(const class_state_infinite &state, int p) {
    Eigen::Tensor<Scalar, 4> temp = get_transfer_matrix_zero(state);
    Eigen::Tensor<Scalar, 4> temp2;
    for(int i = 0; i < p - 2; i++) {
        if(num::mod(i, 2) == 0) {
            temp2 = temp.contract(get_transfer_matrix_LBGA(state), Textra::idx({2, 3}, {0, 1}));

        } else {
            temp2 = temp.contract(get_transfer_matrix_LCGB(state), Textra::idx({2, 3}, {0, 1}));
        }
        temp = temp2;
    }
    return temp;
}

//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_theta(const class_mps_2site  &MPS, Scalar norm)
///*!
// * Returns a two-site MPS
//     @verbatim
//        1--[ LA ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
//                     |                 |
//                     0                 2
//     @endverbatim
// */
//{
//    return
//            MPS.A().contract(MPS.B(), Textra::idx({2},{1})) / norm;
//}
//
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_theta_swapped(const class_mps_2site  &MPS, Scalar norm)
///*!
// * Returns a two-site MPS with A and B swapped
//     @verbatim
//        1--[ LC ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
//                     |                 |
//                     0                 2
//     @endverbatim
// */
//{
//    return MPS.LC_diag() //whatever L_A was in the previous moves
//                    .contract(MPS.B(),        Textra::idx({1},{1}))
//                    .contract(MPS.GA(),       Textra::idx({2},{1}))
//                    .contract(MPS.LA_diag(),       Textra::idx({3},{0}))
//                    .shuffle(Textra::array4{1,0,2,3})
//            /norm;
//}
//
//
//
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_theta_evn(const class_mps_2site  &MPS, Scalar norm)
///*!
// * Returns a right normalized two-site MPS
//     @verbatim
//        1--[ LA ]--[ GA ]-- [ LC ] -- [ GB ]--3
//                     |                 |
//                     0                 2
//     @endverbatim
// */
//
//{
//    return  MPS.A()
//             .contract(MPS.GB(),  Textra::idx({2},{1}))
//            /norm;
//}
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_theta_odd(const class_mps_2site  &MPS, Scalar norm)
///*!
// * Returns a two-site MPS with A and B swapped
//     @verbatim
//        1--[ LC ]--[ GB ]-- [ LB ] -- [ GA ]--3
//                     |                 |
//                     0                 2
//     @endverbatim
// */
//{
//    return MPS.LC_diag()
//                    .contract(MPS.B(),           Textra::idx({1},{1}))
//                    .contract(MPS.GA(),          Textra::idx({2},{1}))
//                    .shuffle(Textra::array4{1,0,2,3})
//            /norm;
//}
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_zero(const class_mps_2site  &MPS) {
//    Eigen::Tensor<Scalar,1> I = MPS.MPS_A->get_LC();
//    I.setConstant(1.0);
//    Eigen::Textra::array<Eigen::IndexPair<long>,0> pair = {};
//
//    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(Textra::array4{0,2,1,3});
//}
//
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_LBGA(const class_mps_2site  &MPS, Scalar norm)  {
//    return MPS.A().contract(MPS.A().conjugate() , Textra::idx({0},{0}))
//                   .shuffle(Textra::array4{0,3,1,2})
//           /norm;
//}
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_GALC(const class_mps_2site  &MPS, Scalar norm)  {
//    return MPS.LC_diag()
//                   .contract(MPS.GA(),               Textra::idx({2},{0}))
//                   .contract(MPS.GA().conjugate(),   Textra::idx({0},{0}))
//                   .contract(MPS.LC_diag(), Textra::idx({3}, {0}) )
//                   .shuffle(Textra::array4{0,2,1,3})
//           /norm;
//}
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_GBLB(const class_mps_2site  &MPS, Scalar norm)  {
//    return MPS.B().contract(MPS.B().conjugate() ,   Textra::idx({0},{0}))
//                   .shuffle(Textra::array4{0,2,1,3})
//           /norm;
//}
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_LCGB(const class_mps_2site  &MPS, Scalar norm)  {
//    return MPS.LC_diag()
//                    .contract(MPS.GB(),               Textra::idx({1},{1}))
//                    .contract(MPS.GB().conjugate(),   Textra::idx({1},{0}))
//                    .contract(MPS.LC_diag(), Textra::idx({2}, {1}) )
//                    .shuffle(Textra::array4{0,3,1,2})
//            /norm;
//}
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_theta_evn(const class_mps_2site  &MPS, Scalar norm)  {
//    using namespace tools::common::views;
//    return get_theta_evn(MPS).contract(get_theta_evn(MPS).conjugate(), Textra::idx({0,2},{0,2})).shuffle(Textra::array4{0,2,1,3}) / norm;
//}
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_theta_odd(const class_mps_2site  &MPS, Scalar norm)  {
//    return get_theta_odd(MPS).contract(get_theta_odd(MPS).conjugate(), Textra::idx({0,2},{0,2})).shuffle(Textra::array4{0,2,1,3}) / norm;
//}
//
//
// Eigen::Tensor<Scalar,4>
// tools::common::views::get_transfer_matrix_AB(const class_mps_2site  &MPS, int p) {
//    Eigen::Tensor<Scalar,4> temp = get_transfer_matrix_zero(MPS);
//    Eigen::Tensor<Scalar,4> temp2;
//    for (int i = 0; i < p-2; i++){
//        if(num::mod(i,2) == 0){
//            temp2 = temp.contract(get_transfer_matrix_LBGA(MPS), Textra::idx({2,3},{0,1}));
//
//        }else{
//            temp2 = temp.contract(get_transfer_matrix_LCGB(MPS), Textra::idx({2,3},{0,1}));
//        }
//        temp = temp2;
//
//
//    }
//    return temp;
//}
//
//
