//
// Created by david on 2019-03-20.
//
//#include <complex.h>
//#undef I

#include "nmspc_quantum_mechanics.h"
#include "general/nmspc_tensor_extra.h"
#include <general/nmspc_iter.h>
#include <io/fmt.h>
#include <io/spdlog.h>
#include <math/linalg.h>
#include <math/num.h>
#include <math/rnd.h>
#include <math/stat.h>
#include <set>
#include <tools/common/log.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

using Scalar = std::complex<double>;
using namespace Eigen;

Eigen::MatrixXcd qm::gen_embedded_spin_operator(const Eigen::MatrixXcd &s, size_t at, size_t sites, bool mirror)
/*
 * Returns a spin operator embedded in a larger Hilbert space. For instance, if at == 1 and sites == 4:
 *
 *   σ¹ = i ⊗ σ ⊗ i ⊗ i
 *
 * where each element is a dxd matrix, resulting in a d^4 * d^4 matrix.

 * Note that if this matrix is converted to a rank-8 tensor, the indexing goes like:
 *
 @verbatim
        3 2 1 0
        | | | |
       [  σ¹  ]
       | | | |
       7 6 5 4
 @endverbatim

 * whereas you would normally want left-to-right indexing in MPS contexts:
 *
 @verbatim
        0 1 2 3
        | | | |
       [  σ¹  ]
       | | | |
       4 5 6 7
 @endverbatim

 * So don't forget to set "reverse = true" if you intend to use the result as a tensor.
 */

{
    if(at >= sites) throw std::logic_error("Expected at < sites. Got [at = " + std::to_string(at) + "] [sites = " + std::to_string(sites) + "]");
    MatrixXcd id     = MatrixXcd::Identity(s.rows(), s.cols());
    MatrixXcd result = at == 0 ? s : id;
    for(size_t site = 1; site < sites; site++)
        result = linalg::matrix::kronecker(result, site == at ? s : id, mirror).eval(); // .eval() is required to avoid aliasing!!
    return result;
}

std::vector<Eigen::MatrixXcd> qm::gen_manybody_spins(const Eigen::MatrixXcd &s, int sites, bool reverse) {
    std::vector<MatrixXcd> S;
    for(int site = 0; site < sites; site++) S.emplace_back(qm::gen_embedded_spin_operator(s, site, sites, reverse));
    return S;
}

namespace qm::spinHalf {

    /* clang-format off */
    Matrix2cd sx = (Matrix2cd() <<
            0.0, 1.0,
            1.0, 0.0).finished();
    Matrix2cd sy = (Matrix2cd() <<
            0.0, imn,
            imp, 0.0).finished();
    Matrix2cd sz = (Matrix2cd() <<
            1.0, 0.0,
            0.0, -1.0).finished();
    Matrix2cd sp = (Matrix2cd() <<
            0.0, 1.0,
            0.0, 0.0).finished();
    Matrix2cd sm = (Matrix2cd() <<
            0.0, 0.0,
            1.0, 0.0).finished();
    Matrix2cd id  = (Matrix2cd() << 1.0, 0.0,
            0.0, 1.0).finished();

    std::array<Vector2cd,2> sx_spinors{(Vector2cd() << 1.0, 1.0).finished()/std::sqrt(2),
                                       (Vector2cd() << 1.0,-1.0).finished()/std::sqrt(2)};

    std::array<Vector2cd,2> sy_spinors{(Vector2cd() << 1.0, imp).finished()/std::sqrt(2),
                                       (Vector2cd() << 1.0, imn).finished()/std::sqrt(2)};

    std::array<Vector2cd,2> sz_spinors{(Vector2cd() << 1.0, 0.0).finished()/std::sqrt(2),
                                       (Vector2cd() << 0.0, 1.0).finished()/std::sqrt(2)};

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
    /* clang-format on */

    Eigen::MatrixXcd gen_embedded_spin_operator(const Eigen::Matrix2cd &s, size_t at, size_t sites, bool swap) {
        // Don't forget to set "swap = true" if you intend to use the result as a tensor.
        return qm::gen_embedded_spin_operator(s, at, sites, swap);
    }

    std::vector<Eigen::Matrix4cd> gen_twobody_spins(const Eigen::Matrix2cd &s, bool swap)
    // Returns a pair of two-body 4x4 spin operators for embedded in a two-site Hilbert space:
    //        (σ ⊗ i, i ⊗ σ)
    // where σ is a 2x2 (pauli) matrix and i is the 2x2 identity matrix.
    // Don't forget to set "swap = true" if you intend to use the result as a tensor.
    {
        std::vector<Matrix4cd> S;
        for(size_t site = 0; site < 2; site++) S.emplace_back(qm::gen_embedded_spin_operator(s, site, 2, swap));
        return S;
    }

}

namespace qm::SpinOne {
    /* clang-format off */

    Matrix3cd sx = (Matrix3cd() <<  0.0, 1.0, 0.0,
                                    1.0, 0.0, 1.0,
                                    0.0, 1.0, 0.0).finished();
    Matrix3cd sy = (Matrix3cd() <<  0.0 , imn, 0.0,
                                    imp,  0.0 ,imn,
                                    0.0 ,  imp, 0.0).finished();
    Matrix3cd sz = (Matrix3cd() <<  1.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0,
                                    0.0, 0.0,-1.0).finished();
    Matrix3cd id = (Matrix3cd()  << 1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0).finished();

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
    /* clang-format on */

    std::vector<Eigen::MatrixXcd> gen_twobody_spins(const Eigen::Matrix3cd &s, bool swap)
    // Returns a pair of two-body 4x4 spin operators for embedded in a two-site Hilbert space:
    //        (σ ⊗ i, i ⊗ σ)
    // where σ is a 3x3 (pauli) matrix and i is the 3x3 identity matrix.
    // So don't forget to set "swap = true" if you intend to use the result as a tensor.
    {
        return qm::gen_manybody_spins(s, 2, swap);
    }
}

namespace qm::timeEvolution {

    std::vector<Eigen::Tensor<Scalar, 2>> Suzuki_Trotter_1st_order(cplx delta_t, const Eigen::Tensor<Scalar, 2> &h_evn, const Eigen::Tensor<Scalar, 2> &h_odd) {
        auto h_evn_matrix = Textra::MatrixMap(h_evn);
        auto h_odd_matrix = Textra::MatrixMap(h_odd);
        return {
            Textra::TensorCast((imn * delta_t * h_evn_matrix).exp()), // exp(-i dt H)
            Textra::TensorCast((imn * delta_t * h_odd_matrix).exp())  // exp(-i dt H)
        };
    }

    std::vector<Eigen::Tensor<Scalar, 2>> Suzuki_Trotter_2nd_order(cplx delta_t, const Eigen::Tensor<Scalar, 2> &h_evn, const Eigen::Tensor<Scalar, 2> &h_odd) {
        auto h_evn_matrix = Textra::MatrixMap(h_evn);
        auto h_odd_matrix = Textra::MatrixMap(h_odd);
        return {Textra::TensorCast((imn * delta_t * h_evn_matrix / 2.0).exp()), Textra::TensorCast((imn * delta_t * h_odd_matrix).exp()),
                Textra::TensorCast((imn * delta_t * h_evn_matrix / 2.0).exp())};
    }

    std::vector<Eigen::Tensor<Scalar, 2>> Suzuki_Trotter_4th_order(cplx delta_t, const Eigen::Tensor<Scalar, 2> &h_evn, const Eigen::Tensor<Scalar, 2> &h_odd)
    /*!
     * Implementation based on
     * Janke, W., & Sauer, T. (1992).
     * Properties of higher-order Trotter formulas.
     * Physics Letters A, 165(3), 199–205.
     * https://doi.org/10.1016/0375-9601(92)90035-K
     *
     */
    {
        auto   h_evn_matrix = Textra::MatrixMap(h_evn);
        auto   h_odd_matrix = Textra::MatrixMap(h_odd);
        double cbrt2        = pow(2.0, 1.0 / 3.0);
        double beta1        = 1.0 / (2.0 - cbrt2);
        double beta2        = -cbrt2 * beta1;
        double alph1        = 0.5 * beta1;
        double alph2        = (1.0 - cbrt2) / 2.0 * beta1;

        std::vector<Eigen::Tensor<Scalar, 2>> temp;
        temp.emplace_back(Textra::TensorCast((alph1 * imn * delta_t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::TensorCast((beta1 * imn * delta_t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::TensorCast((alph2 * imn * delta_t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::TensorCast((beta2 * imn * delta_t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::TensorCast((alph2 * imn * delta_t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::TensorCast((beta1 * imn * delta_t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::TensorCast((alph1 * imn * delta_t * h_evn_matrix).exp()));
        return temp;
    }

    std::vector<Eigen::Tensor<Scalar, 2>> get_twosite_time_evolution_operators(cplx delta_t, size_t susuki_trotter_order, const Eigen::Tensor<Scalar, 2> &h_evn,
                                                                               const Eigen::Tensor<Scalar, 2> &h_odd)
    /*! Returns a set of 2-site unitary gates, using Suzuki Trotter decomposition to order 1, 2 or 3.
     * These gates need to be applied to the MPS one at a time with a swap in between.
     */
    {
        switch(susuki_trotter_order) {
            case 1: return Suzuki_Trotter_1st_order(delta_t, h_evn, h_odd);
            case 2: return Suzuki_Trotter_2nd_order(delta_t, h_evn, h_odd);
            case 4: return Suzuki_Trotter_4th_order(delta_t, h_evn, h_odd);
            default: return Suzuki_Trotter_2nd_order(delta_t, h_evn, h_odd);
        }
    }

    std::vector<Eigen::Tensor<Scalar, 2>> compute_G(const cplx a, size_t susuki_trotter_order, const Eigen::Tensor<Scalar, 2> &h_evn,
                                                    const Eigen::Tensor<Scalar, 2> &h_odd)
    /*! Returns the moment generating function, or characteristic function (if a is imaginary) for the Hamiltonian as a rank 2 tensor.
     *  The legs contain two physical spin indices each
    *   G := exp(iaM) or exp(aM), where a is a small parameter and M is an MPO.
    *   Note that G(-a) = G(a)* if  exp(iaM) !
    *
    @verbatim
                     0
                     |
                [ exp(aH) ]
                     |
                     1
    @endverbatim
    */
    {
        tools::log->warn("compute_G(...): Convention has changed: delta_t, or a, are now multiplied by [-i] in exponentials."
                         " This function may not have been adjusted to the new convention");
        return get_twosite_time_evolution_operators(a, susuki_trotter_order, h_evn, h_odd);
    }

    std::pair<std::vector<qm::Gate>, std::vector<qm::Gate>> get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite) {
        /* Here we do a second-order Suzuki-Trotter decomposition which holds for n-site hamiltonians as described
         * here https://tensornetwork.org/mps/algorithms/timeevo/tebd.html
         * For instance,
         *      H = Sum_a^n H_a
         * where each H_a is a sum of n-site terms.
         *
         * The second-order Suzuki-Trotter decomposition them becomes
         *
         * U2(d) = Prod_{a=1}^n exp(-i[d/2]H_a) Prod_{a=n}^1 exp(-i[d/2]H_a)
         *
         * So this is just the layers applied in reversed order!
         * We return these as a pair of gate layers, and both need to be applied normally for the time evolution
         * to take place
         *
         */

        std::vector<Gate> time_evolution_gates_forward;
        std::vector<Gate> time_evolution_gates_reverse;
        time_evolution_gates_forward.reserve(hams_nsite.size());
        time_evolution_gates_reverse.reserve(hams_nsite.size());

        // Generate first forward layer
        for(auto &h : hams_nsite) {
            time_evolution_gates_forward.emplace_back(h.exp(imn * delta_t * 0.5)); // exp(-i * delta_t * h)
        }
        // Generate second reversed layer
        for(auto &h : iter::reverse(hams_nsite)) {
            time_evolution_gates_reverse.emplace_back(h.exp(imn * delta_t * 0.5)); // exp(-i * delta_t * h)
        }

        // Sanity checks
        if(std::imag(delta_t) == 0) {
            for(auto &t : time_evolution_gates_forward)
                if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                    throw std::runtime_error(fmt::format("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op)));
                }
            for(auto &t : time_evolution_gates_reverse)
                if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                    throw std::runtime_error(fmt::format("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op)));
                }
        }

        return std::make_pair(time_evolution_gates_forward, time_evolution_gates_reverse);
    }

}

std::vector<qm::Gate> qm::lbit::get_unitary_2gate_layer(size_t sites, double fmix) {
    /*! Returns a set of unitary two site operators used to transform between physical and l-bit representations
     *
     *
     @verbatim
                   0      1                0
                   |      |                |
                 [ exp(-ifH) ]  ---> [ exp(-ifH) ]
                   |      |               |
                   2      3               1
     @endverbatim
    */

    tools::log->trace("Generating twosite unitaries");
    constexpr bool kroneckerSwap = false;
    auto           SZ            = qm::spinHalf::gen_twobody_spins(qm::spinHalf::sz, kroneckerSwap); // We use these as matrices
    auto           SP            = qm::spinHalf::gen_twobody_spins(qm::spinHalf::sp, kroneckerSwap); // We use these as matrices
    auto           SM            = qm::spinHalf::gen_twobody_spins(qm::spinHalf::sm, kroneckerSwap); // We use these as matrices
    auto           ID            = qm::spinHalf::gen_twobody_spins(qm::spinHalf::id, kroneckerSwap); // We use these as matrices
    auto           N             = std::vector<Eigen::Matrix4cd>{0.5 * (ID[0] + SZ[0]), 0.5 * (ID[1] + SZ[1])};
    auto           spin_dims     = std::vector<long>{2l,2l};
    std::vector<qm::Gate> unitaries;
    unitaries.reserve(sites - 1);
    for(size_t idx = 0; idx < sites - 1; idx++) {
        double               th0 = rnd::uniform_double_box(1, -1);
        double               th1 = rnd::uniform_double_box(1, -1);
        double               th2 = rnd::uniform_double_box(1, -1);
        double               th3 = rnd::uniform_double_box(1, -1);
        std::complex<double> t(rnd::uniform_double_box(1, -1), rnd::uniform_double_box(1, -1));
        auto indices = std::vector<size_t>{idx, idx + 1};
        Eigen::Matrix4cd H = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
                             SP[0] * SM[1] * t + SP[1] * SM[0] * std::conj(t);

        if constexpr(kroneckerSwap) {
            // Here the kronecker already has index pattern left-to-right and there is no need to shuffle

            //         0               0      1
            //         |               |      |
            //   [ exp(-ifH) ]  ==  [ exp(-ifH) ]
            //        |               |      |
            //        1               2      3

            unitaries.emplace_back(Textra::TensorMap((imn * fmix * H).exp().eval()), indices, spin_dims);
        } else {
            // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
            // the kronecker product that generated two-site gates above has indexed right-to-left
            //         0                   1      0              0      1               0
            //         |                   |      |              |      |               |
            //   [ exp(-ifH) ]  --->    [ exp(-ifH) ]   --->  [ exp(-ifH) ]  --->  [ exp(-ifH) ]
            //        |                   |      |              |      |                |
            //        1                   3      2              2      3                1
            Eigen::Tensor<Scalar, 2> H_shuffled = Textra::TensorMap(H, 2, 2, 2, 2).shuffle(Textra::array4{1, 0, 3, 2}).reshape(Textra::array2{4, 4});
            Eigen::MatrixXcd         expifH     = (imn * fmix * Textra::MatrixMap(H_shuffled)).exp();
            unitaries.emplace_back(Textra::TensorMap(expifH), indices, spin_dims);
        }
    }
    // Sanity check
    for(const auto &u : unitaries)
        if(not Textra::MatrixMap(u.op).isUnitary()) throw std::logic_error("u is not unitary!");
    return unitaries;
}

Eigen::Tensor<Scalar, 2> qm::lbit::get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<Scalar, 2> &H) {
    // Given a matrix H, this returns exp(delta_t * H)
    // For time evolution, just make sure delta_t = -i*d,  where d is a (small) real positive number.
    return Textra::TensorCast((delta_t * Textra::MatrixMap(H)).exp());
}

std::vector<qm::Gate> qm::lbit::get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite) {
    std::vector<Gate> time_evolution_gates;
    time_evolution_gates.reserve(hams_nsite.size());
    for(auto &h : hams_nsite) time_evolution_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
    for(auto &t : time_evolution_gates)
        if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
            throw std::runtime_error(fmt::format("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op)));
        }
    return time_evolution_gates;
}

std::vector<Eigen::Tensor<Scalar, 2>> qm::lbit::get_time_evolution_operators_2site(size_t sites, cplx delta_t,
                                                                                   const std::vector<Eigen::Tensor<Scalar, 2>> &hams_2site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{j,j+1} and h_{j+1,j+2} commute. Therefore we can immediately use the relation
    //      exp(-i*dt *[h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}]) =  exp(-i*dt [h_{j,j+1}]) * exp(-i*dt*[h_{j+1,j+2}]) * ... * exp(-i*dt*[h_{L-2, L-1}])
    // without passing through the Suzuki-Trotter decomposition.
    // Here we expect "hams_2site" to contain terms like  h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}.

    if(hams_2site.size() != sites - 1)
        throw std::logic_error(fmt::format("Wrong number of twosite hamiltonians: {}. Expected {}", hams_2site.size(), sites - 1));

    std::vector<Eigen::Tensor<Scalar, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_2site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

std::vector<Eigen::Tensor<Scalar, 2>> qm::lbit::get_time_evolution_operators_3site(size_t sites, cplx delta_t,
                                                                                   const std::vector<Eigen::Tensor<Scalar, 2>> &hams_3site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{i,j,k} and h_{l,m,n} commute. Therefore we can immediately use the relation
    // exp(A + B) = exp(A)exp(B)
    // without passing through the Suzuki-Trotter decomposition.

    if(hams_3site.size() != sites - 2)
        throw std::logic_error(fmt::format("Wrong number of three-site hamiltonians: {}. Expected {}", hams_3site.size(), sites - 2));

    std::vector<Eigen::Tensor<Scalar, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_3site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

qm::Scalar qm::lbit::get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &tau, size_t pos_tau,
                                        const Eigen::Matrix2cd &sig, size_t pos_sig) {
    // Generate gates for the operators
    tools::log->trace("Computing Trace (tau_{} sig_{})", pos_tau, pos_sig);
//    auto tau_gate = qm::Gate{tau, {pos_tau}, {2l}};
//    auto sig_gate = qm::Gate{sig, {pos_sig}, {2l}};
    auto tau_gate = qm::Gate(tau, {pos_tau}, {2l});
    auto sig_gate = qm::Gate(sig, {pos_sig}, {2l});
    auto                     g   = tau_gate; // Start with the bottom tau gate
    auto                     lc2 = qm::get_lightcone_intersection(unitary_layers, pos_tau, pos_sig);
    bool                     deb = tools::log->level() <= spdlog::level::debug;
    std::vector<std::string> net; // Great for debugging
    std::vector<std::string> log; // Great for debugging
    std::string              empty_layer;
    size_t                   uw = 7; // Width of a unitary 2-site gate box
    size_t                   hw = 3; // Half-width of a unitary 2-site gate box
    size_t                   op = 3; // Overlap of a unitary 2-site gate box
    size_t                   tw = 6; // Tag width
    for(const auto &[idx_layer, layer] : iter::enumerate(unitary_layers)) {
        // Generate
        if(layer.empty()) continue;
        size_t gate_size     = layer.front().pos.size();
        size_t pos_max       = std::max(layer.front().pos.back(), layer.back().pos.back());
        auto   gate_sequence = qm::get_gate_sequence(layer);
        if(deb and net.empty()) {
            auto str_tau = fmt::format("[{}]", pos_tau);
            empty_layer  = fmt::format("{0:^{1}}", " ", tw + pos_max * (uw - op) + op);
            net.emplace_back(empty_layer);
            net.back().replace(0, tw, "tau  :");
            net.back().replace(tw + pos_tau * (uw - op), str_tau.size(), str_tau);
            log.emplace_back(fmt::format("insert tau[{}] now{}", pos_tau, g.pos));
        }
        for(const auto &[idx_sublayer, seq] : iter::enumerate(gate_sequence)) {
            std::string layer_str = empty_layer;
            std::string story_str;
            for(const auto &[idx_seq, pos_gate] : iter::enumerate(seq)) {
                auto &u            = layer.at(pos_gate);
                auto  idx_sublayer = num::mod<size_t>(pos_gate, gate_size);
                if(deb) layer_str.replace(0, tw, fmt::format("u[{:^2}]:", 2 * idx_layer + idx_sublayer)); // Setup layer tag

                // Going through the sequence first forward, then backward, we are handed u gates which may or may not connect to our current g gate.
                // For a successful connection, at least one pos in g should be present in u gate.
                // After connecting a full layer, trace away legs of g that are outside of the light-cone intersection between tau and sigma.

                // Check if g.pos and u.pos have sites in common
                std::vector<size_t> pos_isect;
                std::set_intersection(g.pos.begin(), g.pos.end(), u.pos.begin(), u.pos.end(), back_inserter(pos_isect));

                if(not pos_isect.empty()) {
                    // Found a matching u. Connect it
                    auto pos_old = g.pos;
                    g            = g.insert(u);
                    if(deb)
                        layer_str.replace(tw + u.pos.front() * (uw - op), uw,
                                          fmt::format("[{1:^{0}}]", uw - 2, fmt::format("{:<2},{:>2}", u.pos.front(), u.pos.back())));
                    if(deb) story_str.append(fmt::format("insert u{} ", u.pos));
                }
            }
            // Determine the positions that are allowed
            const std::vector<size_t> &pos_needed = lc2[2 * idx_layer + idx_sublayer + 1]; // This specifies sites that are needed to connect the coming gate
                                                                                           //            if(not pos_needed.empty()){
                                                                                           //
                                                                                           //            }
            // Check if g.pos has non-needed sites
            std::vector<size_t> pos_outside;
            std::set_difference(g.pos.begin(), g.pos.end(), pos_needed.begin(), pos_needed.end(), back_inserter(pos_outside));
            if(not pos_outside.empty()) {
                // Found positions outside of the light cone. Trace them
                auto pos_old = g.pos;
                g            = g.trace_pos(pos_outside);
                for([[maybe_unused]] const auto &p : pos_outside) g.op = g.op * g.op.constant(0.5); // Normalize
                if(deb) story_str.append(fmt::format("trace{} ", pos_outside));
            }
            if(deb) story_str.append(fmt::format("now{} ", g.pos));
            if(deb) net.emplace_back(layer_str);
            if(deb) log.emplace_back(story_str);
        }
    }
    Scalar result;
    if(g.pos.empty()) {
        if(g.op.dimension(0) * g.op.dimension(1) != 1)
            throw std::runtime_error(fmt::format("Expected empty gate to have scalar op: Got dims {}", g.op.dimensions()));
        result = g.op.coeff(0);
        if(deb) net.emplace_back(empty_layer);
        if(deb) net.back().replace(0, tw, "sig  :");
        if(deb) log.emplace_back(fmt::format("sigma not connected -> result = {:.1f}{:+.1f}i", result.real(), result.imag()));
    } else {
        // In the last step we connect the sigma operator and trace everything down to a scalar
        auto num_traces = g.pos.size();
        result          = g.connect_under(sig_gate).trace();
        result *= std::pow(0.5, num_traces); // Normalize
        if(deb) {
            auto str_sig = fmt::format("[{}]", pos_sig);
            net.emplace_back(empty_layer);
            net.back().replace(0, tw, "sig  :");
            net.back().replace(tw + pos_sig * (uw - op), str_sig.size(), str_sig);
            log.emplace_back(fmt::format("insert sigma[{0}] now{1} trace{1} result = {2:.1f}{3:+.1f}i", pos_sig, g.pos, result.real(), result.imag()));
        }
    }

    tools::log->debug("Computed Trace (tau_{} sig_{}) = {:.6f}{:+.6f}i", pos_tau, pos_sig, result.real(), result.imag());
    if(deb)
        for(const auto &[idx, layer] : iter::enumerate_reverse(net)) tools::log->debug("{} | log: {}", layer, log[idx]);

    return result;
}

Eigen::Tensor<qm::Scalar, 2> qm::lbit::get_lbit_real_overlap(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites) {
    Eigen::Tensor<qm::Scalar, 2> lbit_overlap;
    lbit_overlap.resize(static_cast<long>(sites), static_cast<long>(sites));
    for(long j = 0; j < lbit_overlap.dimension(1); j++)
        for(long i = 0; i < lbit_overlap.dimension(0); i++)
            lbit_overlap(i, j) =
                qm::lbit::get_lbit_exp_value(unitary_layers, qm::spinHalf::sz, static_cast<size_t>(i), qm::spinHalf::sz, static_cast<size_t>(j));
    return lbit_overlap;
}

Eigen::Tensor<Scalar, 2> qm::lbit::get_lbit_overlap_permuted(const Eigen::Tensor<Scalar, 2> &lbit_overlap) {
    // First, subtract the center position of each lbit, so we get L lbits centered around zero.
    // In practice, we make a cyclic permutation of the rows of lbit_overlap
    // In addition, we mirror the lbit along its vertical, so that we can average its left and right half together
    long                     rows = lbit_overlap.dimension(0);
    long                     cols = lbit_overlap.dimension(1);
    Eigen::Tensor<Scalar, 2> lbit_overlap_permuted(rows, cols);
    lbit_overlap_permuted.setConstant(0);
    for(long j = 0; j < cols; j++) {
        for(long i = 0; i < rows; i++) {
            // If i < j, this corresponds to the left side of an lbit.
            // Consider i == j to be "origo" for an lbit, so mod(i+j,cols) becomes the index starting from that origo.
            // In addition, we fold the left side of an lbit back onto the right side.
            // To get the correct average, add just half of the value;
            long j_twin   = j;
            long distance = std::abs(i - j);
            long j_perm   = distance;
            if(j >= i) {
                j_twin = i - distance;
                if(j_twin < 0) j_twin = j;
            } else {
                j_twin = i + distance;
                if(j_twin >= cols) j_twin = j;
            }
            tools::log->debug("({:>2},{:>2}): twin ({:>2},{:>2}) | perm = ({:>2},{:>2}) | values {:.6f}, {:.6f} = {:.6f}", i, j, i, j_twin, i, j_perm,
                              std::real(lbit_overlap(i, j)), std::real(lbit_overlap(i, j_twin)),
                              0.5 * (std::real(lbit_overlap(i, j)) + std::real(lbit_overlap(i, j_twin))));
            if(std::abs(std::imag(lbit_overlap(i, j))) > 1e-12)
                tools::log->warn("lbit_overlap({},{}) has imaginary component : |Im({})| > 1e-12", i, j, lbit_overlap(i, j));
            lbit_overlap_permuted(i, j_perm) = 0.5 * (lbit_overlap(i, j) + lbit_overlap(i, j_twin));
        }
    }
    return lbit_overlap_permuted;
}

std::tuple<double, double, std::vector<double>, size_t> qm::lbit::get_characteristic_length_scale(const Eigen::Tensor<Scalar, 2> &lbit_overlap_permuted) {
    // Average along each column to get an estimate of the lbit
    Eigen::Tensor<double, 1> lbit_overlap_avg = lbit_overlap_permuted.real().mean(std::array<long, 1>{0});
    Eigen::Tensor<double, 1> lbit_overlap_log = lbit_overlap_avg.log();
    // Data becomes noisy if the exponential has decayed, so find a cutoff to get the slope using only the first part of the curve
    auto y            = std::vector<double>(lbit_overlap_log.data(), lbit_overlap_log.data() + lbit_overlap_log.size());
    auto x            = num::range<double>(0, y.size());
    auto v            = stat::find_last_valid_point(y);
    auto c            = std::count_if(y.begin(), y.begin() + v, [](auto &val) { return val > -12.0; });
    auto [slope, res] = stat::slope(x, y, 0, c);
    //    tools::log->debug("lbit overlap : \n{}", linalg::tensor::to_string(lbit_overlap.real(), 6));
    //    tools::log->debug("lbit permuted: \n{}", linalg::tensor::to_string(lbit_overlap_permuted, 6));
    //    tools::log->debug("lbit averaged: \n{}", linalg::tensor::to_string(lbit_overlap_average, 6));
    //    tools::log->debug("lbit logged  : \n{}", linalg::tensor::to_string(lbit_overlap_log, 6));
    //    tools::log->debug("lbit vectored: \n{}", y);
    double cls = 1.0 / std::abs(slope);
    tools::log->debug("Computed lbit width {:.6f} | sse {:.6f} | using {} points", cls, res, c);
    //    auto yavg = std::vector<double>(lbit_overlap_avg.data(),lbit_overlap_avg.data() + c);
    auto yavg = std::vector<double>(lbit_overlap_avg.data(), lbit_overlap_avg.data() + lbit_overlap_avg.size());
    return {cls, res, yavg, c};
}

Eigen::Tensor<Scalar, 2> qm::lbit::get_lbit_overlap_averaged(const std::vector<Eigen::Tensor<Scalar, 2>> &lbit_overlap_vec) {
    Eigen::Tensor<Scalar, 2> avg;
    //    Eigen::Tensor<Scalar,2> err;
    if(not lbit_overlap_vec.empty()) {
        long                rows = lbit_overlap_vec.front().dimension(0);
        long                cols = lbit_overlap_vec.front().dimension(1);
        size_t              reps = lbit_overlap_vec.size();
        std::vector<Scalar> slice(reps);
        avg.resize(rows, cols);
        for(long c = 0; c < cols; c++) {
            for(long r = 0; r < rows; r++) {
                for(const auto &[i, elem] : iter::enumerate(lbit_overlap_vec)) {
                    if(std::abs(std::imag(elem(r, c))) > 1e-12)
                        tools::log->warn("elem[{}]({},{}) has imaginary component : |Im({})| > 1e-12", i, r, c, elem(r, c));
                    slice[i] = elem(r, c);
                }
                avg(r, c) = stat::mean(slice);
                //                err(r,c) = stat::sterr(slice);
            }
        }
    }
    return avg;
    //    return {avg,err};
}
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::Tensor<double, 3>,Eigen::Tensor<double, 4>>
    qm::lbit::get_lbit_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec, size_t sites, size_t reps) {
    long                     rows = static_cast<long>(fmix_vec.size());
    long                     cols = static_cast<long>(udepth_vec.size());
    Eigen::MatrixXd          cls_avg(rows, cols);
    Eigen::MatrixXd          cls_err(rows, cols);
    Eigen::MatrixXd          sse_avg(rows, cols);
    Eigen::MatrixXd          sse_err(rows, cols);
    Eigen::Tensor<double, 3> lbit_decay(rows, cols, sites);
    Eigen::Tensor<double, 4> lbit_lioms(rows, cols, sites, sites);
    lbit_decay.setZero();
    lbit_lioms.setZero();
    std::array<long, 3> offset3{}, extent3{};
    std::array<long, 4> offset4{}, extent4{};

#pragma omp parallel for collapse(2) schedule(dynamic)
    for(size_t uidx = 0; uidx < udepth_vec.size(); uidx++) {
        for(size_t fidx = 0; fidx < fmix_vec.size(); fidx++) {
            //    for(const auto & [uidx,udep] : iter::enumerate(udepth_vec)){

            //        for (const auto & [fidx,fmix] : iter::enumerate(fmix_vec) ){
            auto                                  fmix = fmix_vec[fidx];
            auto                                  udep = udepth_vec[uidx];
            std::vector<double>                   cls_vec(reps);
            std::vector<double>                   sse_vec(reps);
            std::vector<Eigen::Tensor<Scalar, 2>> lbit_overlap_vec(reps);
            for(size_t i = 0; i < reps; i++) {
                std::vector<std::vector<qm::Gate>> layers;
                for(size_t l = 0; l < udep; l++) layers.emplace_back(qm::lbit::get_unitary_2gate_layer(sites, fmix));
                lbit_overlap_vec[i] = qm::lbit::get_lbit_real_overlap(layers, sites);
            }

            auto lbit_overlap_avg = qm::lbit::get_lbit_overlap_averaged(lbit_overlap_vec);
            auto lbit_overlap_per = qm::lbit::get_lbit_overlap_permuted(lbit_overlap_avg);

            auto [cls, sse, y, c] = qm::lbit::get_characteristic_length_scale(lbit_overlap_per);
#if defined(_OPENMP)
            tools::log->info("Computed u {} | f {:.4f} | lbit width {:.6f} | sse {:.6f} | threads {} | points {}: {:.8f}", udep, fmix, cls, sse,
                             omp_get_num_threads(), c, fmt::join(y, ", "));
#else
            tools::log->info("Computed u {} | f {:.4f} | lbit width {:.6f} | sse {:.6f} | points {}: {:.8f}", udep, fmix, cls, sse, c, fmt::join(y, ", "));
#endif
            cls_avg(static_cast<long>(fidx), static_cast<long>(uidx)) = cls;
            sse_avg(static_cast<long>(fidx), static_cast<long>(uidx)) = sse;

            offset3                             = {static_cast<long>(fidx), static_cast<long>(uidx), 0};
            extent3                             = {1, 1, static_cast<long>(y.size())};
            offset4                             = {static_cast<long>(fidx), static_cast<long>(uidx), 0, 0};
            extent4                             = {1, 1, lbit_overlap_avg.dimension(0), lbit_overlap_avg.dimension(1)};
            lbit_decay.slice(offset3, extent3) = Eigen::TensorMap<Eigen::Tensor<double, 1>>(y.data(), y.size());
            lbit_lioms.slice(offset4, extent4) = lbit_overlap_avg.real().reshape(extent4);
        }
    }
    return {cls_avg, sse_avg, lbit_decay, lbit_lioms};
}

std::tuple<Eigen::Tensor<Scalar, 4>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>> qm::mpo::pauli_mpo(const Eigen::MatrixXcd &paulimatrix)
/*! Builds the MPO string for measuring  spin on many-body systems.
 *      P = Π  s_{i}
 * where Π is the product over all sites, and s_{i} is the given pauli matrix for site i.
 *
 * MPO = | s | (a 1 by 1 matrix with a single pauli matrix element)
 *
 *        2
 *        |
 *    0---s---1
 *        |
 *        3
 *
 */
{
    long                     spin_dim = paulimatrix.rows();
    std::array<long, 4>      extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>      extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO(1, 1, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrix);

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 1); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 1); // The right edge
    Ledge(0, 0, 0) = 1;
    Redge(0, 0, 0) = 1;
    return std::make_tuple(MPO, Ledge, Redge);
}

std::tuple<std::vector<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::parity_projector_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites, int sign)
/*! Builds the MPO that projects out the MPS component in a parity sector.
 * |psi+->  = P |psi>=  1/2 (1 +- S) |psi>
 * Here 1 = outer product of L=sites 2x2 identity matrices, i.e. Kron_(i=0)^(L-1) I_(2x2)
 * Also S = outer product of L=sites 2x2 pauli matrices, i.e. Kron_(i=0)^(L-1) s_(2x2)
 * The sign and the factor 1/2 is put into the left edge at the end.
 *
 *                     | I   0  |
 *    S   =      1/2 * | 0   s  |
 *
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    long                     spin_dim = paulimatrix.rows();
    auto                     I        = Eigen::MatrixXcd::Identity(spin_dim, spin_dim).eval();
    std::array<long, 4>      extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>      extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO(2, 2, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(I);
    MPO.slice(std::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrix);

    std::vector<Eigen::Tensor<Scalar, 4>> mpos(sites, MPO);
    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 2); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 2); // The right edge
    Ledge(0, 0, 0) = 0.5;                    // 0.5;
    Ledge(0, 0, 1) = 0.5 * sign;
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;

    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::vector<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i, where S and I is a pauli matrix or an identity matrix respectively
 *
 * MPO = | s |
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    long                     spin_dim = paulimatrix.rows();
    std::array<long, 4>      extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>      extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_I(1, 1, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_S(1, 1, spin_dim, spin_dim);
    MPO_I.setZero();
    MPO_S.setZero();
    MPO_I.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrix);
    MPO_S.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorCast(Eigen::MatrixXcd::Identity(spin_dim, spin_dim));

    // We have to push in an even number of pauli matrices to retain the parity sector.
    // Choosing randomly
    std::vector<int> binary(sites, -1);
    int              sum = 0;
    while(true) {
        binary[rnd::uniform_integer_box<size_t>(0, sites)] *= -1;
        sum = std::accumulate(binary.begin(), binary.end(), 0);
        if((num::mod<size_t>(sites, 2) == 0 and sum == 0) or (num::mod<size_t>(sites, 2) == 1 and sum == 1)) break;
    }

    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    for(auto &val : binary) {
        if(val < 0)
            mpos.push_back(MPO_S);
        else
            mpos.push_back(MPO_I);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 1); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 1); // The right edge
    Ledge(0, 0, 0) = 1;
    Redge(0, 0, 0) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::vector<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos_x2(const Eigen::MatrixXcd &paulimatrix1, const Eigen::MatrixXcd &paulimatrix2, const size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i.
 * S is the sum of pauli matrices s1 and s2, and where I is an identity matrix of the same size
 *            | s1  0  |
 * S   =      | 0   s2 |
 *
 *            | id  0  |
 * I   =      | 0   id |
 *
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    if(paulimatrix1.rows() != paulimatrix2.rows()) throw std::logic_error("Pauli matrices must be of equal size");
    long                     spin_dim = paulimatrix1.rows();
    auto                     I        = Eigen::MatrixXcd::Identity(spin_dim, spin_dim).eval();
    std::array<long, 4>      extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>      extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_S(2, 2, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_I(2, 2, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_P(2, 2, spin_dim, spin_dim);
    MPO_S.setZero();
    MPO_I.setZero();
    MPO_P.setZero();

    MPO_S.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrix1);
    MPO_S.slice(std::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrix2);
    MPO_I.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(I);
    MPO_I.slice(std::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(I);
    MPO_P.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(I);
    MPO_P.slice(std::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrix1);

    // Push in an even number of operators
    std::vector<int> binary(sites, -1);
    int              sum = 0;
    while(true) {
        binary[rnd::uniform_integer_box<size_t>(0, sites)] *= -1;
        sum = std::accumulate(binary.begin(), binary.end(), 0);
        if((num::mod<size_t>(sites, 2) == 0 and sum == 0) or (num::mod<size_t>(sites, 2) == 1 and sum == 1)) break;
    }
    if(binary.size() != sites) throw std::logic_error("Size mismatch");
    // Generate the list
    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    for(auto &val : binary) {
        if(val < 0)
            mpos.push_back(MPO_S);
        else
            mpos.push_back(MPO_I);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 2); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 2); // The right edge
    Ledge(0, 0, 0) = 1.0 / std::sqrt(2);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(2);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::vector<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i.
 * S is the sum of pauli matrices s0,s1,s2... , and where I is an identity matrix of the same size
 *
 *            | s0  0   0  .  |
 * S   =      | 0   s1  0  .  |
 *            | 0   0  s2  .  |
 *            | .   .   . ... |
 *
 *            | id  0   0  .  |
 * I   =      | 0   id  0  .  |
 *            | 0   0  id  .  |
 *            | .   .   . ... |
 *
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    if(paulimatrices.empty()) throw std::runtime_error("List of pauli matrices is empty");
    long                     num_paulis = static_cast<long>(paulimatrices.size());
    long                     spin_dim   = 2;
    auto                     I          = Eigen::MatrixXcd::Identity(spin_dim, spin_dim).eval();
    std::array<long, 4>      extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>      extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_I(num_paulis, num_paulis, spin_dim, spin_dim);
    MPO_S.setZero();
    MPO_I.setZero();
    for(long diag_pos = 0; diag_pos < num_paulis; diag_pos++) {
        MPO_S.slice(std::array<long, 4>{diag_pos, diag_pos, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(paulimatrices[static_cast<size_t>(diag_pos)]);
        MPO_I.slice(std::array<long, 4>{diag_pos, diag_pos, 0, 0}, extent4).reshape(extent2) = Textra::TensorMap(I);
    }

    // Push in an even number of operators
    // This is so that we get a 50% chance of applying a gate.
    std::vector<int> binary(sites, -1);
    int              sum = 0;
    while(true) {
        binary[rnd::uniform_integer_box<size_t>(0, sites - 1)] *= -1;
        sum = std::accumulate(binary.begin(), binary.end(), 0);
        if((num::mod<size_t>(sites, 2) == 0 and sum == 0) or (num::mod<size_t>(sites, 2) == 1 and sum == 1)) break;
    }
    if(binary.size() != sites) throw std::logic_error("Size mismatch");
    // Generate the list
    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    std::vector<std::string>              mpos_str;
    for(auto &val : binary) {
        if(val < 0) {
            mpos.push_back(MPO_S);
            mpos_str.emplace_back("S");
        } else {
            mpos.push_back(MPO_I);
            mpos_str.emplace_back("I");
        }
    }
    tools::log->warn("Generated random pauli MPO string: {}", mpos_str);
    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, num_paulis); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, num_paulis); // The right edge
    Ledge(0, 0, 0) = 1.0 / std::sqrt(num_paulis);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(num_paulis);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::vector<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::sum_of_pauli_mpo(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites, RandomizerMode mode)
/*! Builds a string of MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i are MPOs with 2x2 (pauli) matrices on the diagonal
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 * If mode == RandomizerMode::SHUFFLE:
 *
 *                 | s0  0   0  .  |
 *      O_i =      | 0   s1  0  .  |
 *                 | 0   0  s2  .  |
 *                 | .   .   . ... |
 *
 * where for each O_i the matrices s0, s1, s2 are shuffled randomly
 *
 * If mode == RandomizerMode::SELECT1:
 *
 *      O_i =  | s  |
 *
 *  where for each O_i one of the matrices s0, s1, s2... is selected randomly
 *
 * If mode == RandomizerMode::ASIS:
 *
 *                 | s0  0   0  .  |
 *      O_i =      | 0   s1  0  .  |
 *                 | 0   0  s2  .  |
 *                 | .   .   . ... |
 *
 * where for each O_i the matrices s0, s1, s2... are placed in order as given
 *
 */

{
    if(paulimatrices.empty()) throw std::runtime_error("List of pauli matrices is empty");
    long                spin_dim = 2;
    std::array<long, 4> extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2> extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    std::array<long, 4> offset4  = {0, 0, 0, 0};

    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    auto                                  pauli_idx = num::range<size_t>(0, paulimatrices.size(), 1);

    for(size_t site = 0; site < sites; site++) {
        Eigen::Tensor<Scalar, 4> mpo;
        switch(mode) {
            case RandomizerMode::SELECT1: {
                mpo.resize(1, 1, spin_dim, spin_dim);
                mpo.setZero();
                auto        rnd_idx                          = rnd::uniform_integer_box<size_t>(0, paulimatrices.size() - 1);
                const auto &pauli                            = paulimatrices[rnd_idx];
                mpo.slice(offset4, extent4).reshape(extent2) = Textra::TensorCast(pauli);
                break;
            }
            case RandomizerMode::SHUFFLE: {
                rnd::shuffle(pauli_idx);
                [[fallthrough]];
            }
            case RandomizerMode::ASIS: {
                auto num_paulis = static_cast<long>(paulimatrices.size());
                mpo.resize(num_paulis, num_paulis, spin_dim, spin_dim);
                for(long idx = 0; idx < num_paulis; idx++) {
                    auto        uidx                             = static_cast<size_t>(idx);
                    const auto &pauli                            = paulimatrices[pauli_idx[uidx]];
                    offset4                                      = {idx, idx, 0, 0};
                    mpo.slice(offset4, extent4).reshape(extent2) = Textra::TensorCast(pauli);
                }
                break;
            }
        }
        mpos.emplace_back(mpo);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 1); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 1); // The right edge
    switch(mode) {
        case RandomizerMode::SHUFFLE:
        case RandomizerMode::ASIS: {
            Ledge.resize(1, 1, paulimatrices.size());
            Redge.resize(1, 1, paulimatrices.size());
            for(size_t idx = 0; idx < paulimatrices.size(); idx++) {
                Ledge(0, 0, idx) = 1.0 / std::sqrt(paulimatrices.size());
                Redge(0, 0, idx) = 1;
            }
            break;
        }
        case RandomizerMode::SELECT1: {
            Ledge(0, 0, 0) = 1;
            Redge(0, 0, 0) = 1;
            break;
        }
    }
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::vector<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, const std::vector<double> &uniform_dist_widths, size_t sites)
/*! Builds a set of MPO's used for randomizing a state  pauli matrix MPO's with random weights picked from a uniform distribution
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is the MPO sum of pauli matrices with random weights.
 *
 *            | c0*s0   0       0     .   |
 * O_i =      | 0       c1*s1   0     .   |
 *            | 0       0       c2*s2 .   |
 *            | .       .       .     ... |
 *  Here s_i are 2x2 pauli matrices (including identity) and
 *  the weight coefficients c_i are random real numbers drawn from a uniform distribution U(-w,w).
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    if(paulimatrices.empty()) throw std::runtime_error("List of pauli matrices is empty");
    if(paulimatrices.size() != uniform_dist_widths.size()) throw std::runtime_error("List size mismatch: paulimatrices and uniform_dist_widths");
    long                num_paulis = static_cast<long>(paulimatrices.size());
    long                spin_dim   = 2;
    std::array<long, 4> extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2> extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */

    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    for(size_t site = 0; site < sites; site++) {
        Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
        MPO_S.setZero();
        for(long idx = 0; idx < num_paulis; idx++) {
            auto        uidx                               = static_cast<size_t>(idx);
            auto        coeff                              = 1 + rnd::uniform_double_box(uniform_dist_widths[uidx]);
            auto        offset4                            = std::array<long, 4>{idx, idx, 0, 0};
            const auto &pauli                              = paulimatrices[uidx];
            MPO_S.slice(offset4, extent4).reshape(extent2) = Textra::TensorCast(coeff * pauli);
        }
        mpos.emplace_back(MPO_S);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, num_paulis); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, num_paulis); // The right edge
    Ledge(0, 0, 0) = 1.0 / std::sqrt(num_paulis);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(num_paulis);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

//[ 8   0   0   0   0   0   0   0 ]
//[ 0  16   0   0   0   0   0   0 ]
//[ 0   0  32   0   0   0   0   0 ]
//[ 0   0   0  32   0   0   0   0 ]
//[ 0   0   0   0  32   0   0   0 ]
//[ 0   0   0   0   0  32   0   0 ]
//[ 0   0   0   0   0   0  16   0 ]
//[ 0   0   0   0   0   0   0   8 ]