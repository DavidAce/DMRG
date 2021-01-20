//
// Created by david on 2019-03-20.
//
//#include <complex.h>
//#undef I

#include "nmspc_quantum_mechanics.h"
#include "general/nmspc_tensor_extra.h"
#include <Eigen/Core>
#include <general/nmspc_iter.h>
#include <math/num.h>
#include <math/rnd.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <set>
#include <vector>

using Scalar = std::complex<double>;
using namespace Eigen;

Eigen::MatrixXcd qm::gen_embedded_spin_operator(const Eigen::MatrixXcd &s, size_t at, size_t sites, bool swap)
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

 * whereas you would normally want left-to-right indexing:
 *
 @verbatim
        0 1 2 3
        | | | |
       [  σ¹  ]
       | | | |
       4 5 6 7
 @endverbatim

 * So don't forget to set "swap = true" if you intend to use the result as a tensor.
 */

{
    if(at >= sites) throw std::logic_error("Expected at < sites. Got [at = " + std::to_string(at) + "] [sites = " + std::to_string(sites) + "]");
    MatrixXcd id     = MatrixXcd::Identity(s.rows(), s.cols());
    MatrixXcd result = at == 0 ? s : id;
    if(swap)
        for(size_t site = 1; site < sites; site++) result = kroneckerProduct(site == at ? s : id, result).eval(); // .eval() is required to avoid aliasing!!
    else
        for(size_t site = 1; site < sites; site++) result = kroneckerProduct(result, site == at ? s : id).eval(); // .eval() is required to avoid aliasing!!
    return result;
}

std::vector<Eigen::MatrixXcd> qm::gen_manybody_spins(const Eigen::MatrixXcd &s, int sites, bool swap) {
    std::vector<MatrixXcd> S;
    for(int site = 0; site < sites; site++) S.emplace_back(qm::gen_embedded_spin_operator(s, site, sites, swap));
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
        auto h_evn_matrix = Textra::TensorMatrixMap(h_evn);
        auto h_odd_matrix = Textra::TensorMatrixMap(h_odd);
        return {
            Textra::MatrixToTensor((imn * delta_t * h_evn_matrix).exp()), // exp(-i dt H)
            Textra::MatrixToTensor((imn * delta_t * h_odd_matrix).exp())  // exp(-i dt H)
        };
    }

    std::vector<Eigen::Tensor<Scalar, 2>> Suzuki_Trotter_2nd_order(cplx delta_t, const Eigen::Tensor<Scalar, 2> &h_evn, const Eigen::Tensor<Scalar, 2> &h_odd) {
        auto h_evn_matrix = Textra::TensorMatrixMap(h_evn);
        auto h_odd_matrix = Textra::TensorMatrixMap(h_odd);
        return {Textra::MatrixToTensor((imn * delta_t * h_evn_matrix / 2.0).exp()), Textra::MatrixToTensor((imn * delta_t * h_odd_matrix).exp()),
                Textra::MatrixToTensor((imn * delta_t * h_evn_matrix / 2.0).exp())};
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
        auto   h_evn_matrix = Textra::TensorMatrixMap(h_evn);
        auto   h_odd_matrix = Textra::TensorMatrixMap(h_odd);
        double cbrt2        = pow(2.0, 1.0 / 3.0);
        double beta1        = 1.0 / (2.0 - cbrt2);
        double beta2        = -cbrt2 * beta1;
        double alph1        = 0.5 * beta1;
        double alph2        = (1.0 - cbrt2) / 2.0 * beta1;

        std::vector<Eigen::Tensor<Scalar, 2>> temp;
        temp.emplace_back(Textra::MatrixToTensor((alph1 * imn * delta_t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((beta1 * imn * delta_t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((alph2 * imn * delta_t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((beta2 * imn * delta_t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((alph2 * imn * delta_t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((beta1 * imn * delta_t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((alph1 * imn * delta_t * h_evn_matrix).exp()));
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

    std::vector<qm::Gate> unitaries;
    unitaries.reserve(sites - 1);
    for(size_t idx = 0; idx < sites - 1; idx++) {
        double               th0 = rnd::uniform_double_box(1,-1);
        double               th1 = rnd::uniform_double_box(1,-1);
        double               th2 = rnd::uniform_double_box(1,-1);
        double               th3 = rnd::uniform_double_box(1,-1);
        std::complex<double> t(rnd::uniform_double_box(1,-1), rnd::uniform_double_box(1,-1));

        Eigen::Matrix4cd     H =
            th3 * N[0] * N[1] +
            th2 * N[1] * (ID[0] - N[0]) +
            th1 * N[0] * (ID[1] - N[1]) +
            th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
            SP[0] * SM[1] * t +
            SP[1] * SM[0] * std::conj(t);

        if constexpr(kroneckerSwap) {
            // Here the kronecker already has index pattern left-to-right and there is no need to shuffle

            //         0               0      1
            //         |               |      |
            //   [ exp(-ifH) ]  ==  [ exp(-ifH) ]
            //        |               |      |
            //        1               2      3

            unitaries.emplace_back(Textra::MatrixToTensor2((imn * fmix * H).exp()), std::vector<size_t>{idx, idx + 1}, std::vector<long>{2,2});
        } else {
            // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
            // the kronecker product that generated two-site gates above has indexed right-to-left
            //         0                   1      0              0      1               0
            //         |                   |      |              |      |               |
            //   [ exp(-ifH) ]  --->    [ exp(-ifH) ]   --->  [ exp(-ifH) ]  --->  [ exp(-ifH) ]
            //        |                   |      |              |      |                |
            //        1                   3      2              2      3                1
            Eigen::Tensor<Scalar, 2> H_shuffled = Textra::MatrixTensorMap(H, 2, 2, 2, 2).shuffle(Textra::array4{1, 0, 3, 2}).reshape(Textra::array2{4, 4});
            Eigen::MatrixXcd         expifH     = (imn * fmix * Textra::TensorMatrixMap(H_shuffled)).exp();
            unitaries.emplace_back(Textra::MatrixTensorMap(expifH), std::vector<size_t>{idx, idx + 1}, std::vector<long>{2, 2});
        }
    }
    // Sanity check
    for(const auto &u : unitaries)
        if(not Textra::TensorMatrixMap(u.op).isUnitary()) throw std::logic_error("u is not unitary!");
    return unitaries;
}

Eigen::Tensor<Scalar, 2> qm::lbit::get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<Scalar, 2> &H) {
    // Given a matrix H, this returns exp(delta_t * H)
    // For time evolution, just make sure delta_t = -i*d,  where d is a (small) real positive number.
    return Textra::MatrixToTensor((delta_t * Textra::TensorMatrixMap(H)).exp());
}

std::vector<qm::Gate> qm::lbit::get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite) {
    std::vector<Gate> time_evolution_gates;
    time_evolution_gates.reserve(hams_nsite.size());
    for(auto &h : hams_nsite) time_evolution_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
    for(auto &t : time_evolution_gates)
        if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
            std::cout << "Supposedly not unitary: \n " << t.op << std::endl;
            throw std::runtime_error(fmt::format("Time evolution operator at pos {} is not unitary", t.pos));
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


std::vector<const qm::Gate*> get_adjacent_u(const std::vector<qm::Gate> & u_layer, const std::vector<size_t> & pos){
   std::vector<const qm::Gate*> adjs;
   adjs.reserve(pos.size());
    for (const auto & u : u_layer){
        if(u.pos.front() > pos.back()) continue; // Skip pos too far
        if(u.pos.back() < pos.front()) continue; // Skip pos too far
        if(u.pos == pos) continue; // Do not include self
        if(std::any_of(u.pos.begin(), u.pos.end(),[& pos](const auto & p){ return std::find(pos.begin(),pos.end(),p) != pos.end();})) adjs.push_back(&u);
        if(adjs.size() >= pos.size()) break;
   }
   return adjs;
}



std::vector<size_t> get_pos_in_lightcone(size_t gate_size, size_t pos_max,size_t num_layers,  size_t idx_layer, size_t idx_sublayer, size_t pos_tau, size_t pos_sig){
    // Note that idx_layer2 and num_layers2 must follow the doubled layer convention.
    // Let 1 normal layer be the set of gates such that each position on the chain appears exactly once on a left-most gate leg.
    // In such a layer, every position appears "gate_size" times in total (each time in a different gate).
    // In this convention, each layer therefore contains "gate_size" number of sublayers.
    // Therefore idx_layer indexes the normal "outer" layer, and idx_sublayer indexes the inner layer.
    long gate_sizel = static_cast<long>(gate_size);
    long step_size      = gate_sizel - 1;
    long step_size_odd  = num::mod(step_size, 2l);
    long idx_layer2_tau = static_cast<long>(2*idx_layer + idx_sublayer);  // Index sublayers counting from tau
    long idx_layer2_sig = static_cast<long>(2*num_layers - (2*idx_layer + idx_sublayer) - 1);  // Index sublayers counting from sig
    long posl_max = static_cast<long>(pos_max);
    long posl_tau = static_cast<long>(pos_tau);
    long posl_sig = static_cast<long>(pos_sig);
    long pos_tau_left = std::clamp(posl_tau - num::mod(posl_tau, gate_sizel), 0l, posl_max); // Get the left-most position of the gate on which tau is attached
    long pos_sig_left = std::clamp(posl_sig - num::mod(posl_sig + step_size_odd, gate_sizel), 0l, posl_max); // Get the left-most position of the gate on which sig is attached
//    long keep_legs = idx_layer2_sig == 0 ? 1 : 0;
//    long sig_offset = idx_layer2_sig == 0 ? 0 : 1;

    long lightcone_pos_tau_min = std::max(static_cast<long>(pos_tau_left) - static_cast<long>(idx_layer2_tau), 0l);
    long lightcone_pos_tau_max = std::min(static_cast<long>(pos_tau_left) + static_cast<long>(idx_layer2_tau + step_size), static_cast<long>(pos_max));
    long lightcone_pos_sig_min = std::max(static_cast<long>(pos_sig_left) - static_cast<long>(idx_layer2_sig), 0l);
    long lightcone_pos_sig_max = std::min(static_cast<long>(pos_sig_left + step_size) + static_cast<long>(idx_layer2_sig), static_cast<long>(pos_max));
    long lightcone_pos_min = std::max(lightcone_pos_tau_min, lightcone_pos_sig_min);
    long lightcone_pos_max = std::min(lightcone_pos_tau_max, lightcone_pos_sig_max);
    if(lightcone_pos_min < 0) throw std::logic_error("lightcone_pos_min < 0");
    if(lightcone_pos_max > static_cast<long>(pos_max)) throw std::logic_error("lightcone_pos_max > pos_max");
    tools::log->info("Lightcone layer [{},{}] = {} ({})| tau [{} {}] sig [{} {}] | lim [{} {}]", idx_layer, idx_sublayer, idx_layer2_tau,idx_layer2_sig, lightcone_pos_tau_min,lightcone_pos_tau_max,lightcone_pos_sig_min,lightcone_pos_sig_max,lightcone_pos_min,lightcone_pos_max  );
    std::vector<size_t> pos;
    long p = lightcone_pos_min;
    while(p <= lightcone_pos_max) pos.emplace_back(static_cast<size_t>(p++));
    if(not pos.empty() and pos.front() != static_cast<size_t>(lightcone_pos_min)) throw std::logic_error("pos.front() != lightcone_pos_min");
    if(not pos.empty() and pos.back()  != static_cast<size_t>(lightcone_pos_max)) throw std::logic_error("pos.back() != lightcone_pos_max");
    return pos;
}


std::vector<std::vector<size_t>> get_pos_in_lightcone2(const std::vector<std::vector<qm::Gate>> & unitary_layers, size_t pos_tau, size_t pos_sig){
    std::vector<std::vector<size_t>> tau_cone, sig_cone, res_cone;
    tau_cone.emplace_back(std::vector<size_t>{pos_tau});
    sig_cone.emplace_back(std::vector<size_t>{pos_sig});
    for(const auto & [idx_layer,layer] : iter::enumerate(unitary_layers)) {
        // Generate
        if(layer.empty()) continue;
        std::vector<std::vector<size_t>> gate_sequence;
        size_t gate_size = layer.front().pos.size();
        size_t pos_max   = layer.back().pos.back();
        for(size_t offset = 0; offset < gate_size; offset++) {
            if(offset + gate_size > pos_max +1) break;
            auto off_idx = num::range<size_t>(offset, pos_max - gate_size + 2, gate_size);
            if(num::mod<size_t>(offset, 2) == 1) std::reverse(off_idx.begin(), off_idx.end()); // If odd, reverse the sequence
            gate_sequence.emplace_back(off_idx);
        }
        for(const auto &[idx_sublayer, seq] : iter::enumerate(gate_sequence)) {
            std::set<size_t> match;
            for(const auto & [idx_seq, pos_gate] : iter::enumerate(seq)){
                auto &u = layer[pos_gate];
                std::vector<size_t> pos_isect;
                std::set_intersection(tau_cone.back().begin(),tau_cone.back().end(),
                                      u.pos.begin(),u.pos.end(),
                                      back_inserter(pos_isect));
                if(not pos_isect.empty())
                    match.insert(u.pos.begin(),u.pos.end()); // Add positions from u if the gate connects to the current cone
            }
            if(match.empty()) tau_cone.emplace_back(tau_cone.back());
            else tau_cone.emplace_back(std::vector<size_t>(match.begin(),match.end()));
        }
    }
    // Do the same in reverse for sig
    for(const auto & [idx_layer,layer] : iter::enumerate_reverse(unitary_layers)) {
        // Generate
        if(layer.empty()) continue;
        std::vector<std::vector<size_t>> gate_sequence;
        size_t gate_size = layer.front().pos.size();
        size_t pos_max   = layer.back().pos.back();
        for(size_t offset = 0; offset < gate_size; offset++) {
            if(offset + gate_size > pos_max +1) break;
            auto off_idx = num::range<size_t>(offset, pos_max - gate_size + 2, gate_size);
            if(num::mod<size_t>(offset, 2) == 1) std::reverse(off_idx.begin(), off_idx.end()); // If odd, reverse the sequence
            gate_sequence.emplace_back(off_idx);
        }
        for(const auto &[idx_sublayer, seq] : iter::enumerate_reverse(gate_sequence)) {
            std::set<size_t> match;
            for(const auto & [idx_seq, pos_gate] : iter::enumerate_reverse(seq)){
                auto &u = layer[pos_gate];
                std::vector<size_t> pos_isect;
                std::set_intersection(sig_cone.back().begin(),sig_cone.back().end(),
                                      u.pos.begin(),u.pos.end(),
                                      back_inserter(pos_isect));
                if(not pos_isect.empty())
                    match.insert(u.pos.begin(),u.pos.end()); // Add positions from u if the gate connects to the current cone
            }
            if(match.empty()) sig_cone.emplace_back(sig_cone.back());
            else sig_cone.emplace_back(std::vector<size_t>(match.begin(),match.end()));
        }
    }
    if(tau_cone.size() != sig_cone.size()) throw std::runtime_error("tau and sig cones should have equal size!");
    // Now turn the sig_cone upside down
    std::reverse(sig_cone.begin(),sig_cone.end());
    // Find the intersection between tau and sig cones
    for(size_t idx_sublayer = 0; idx_sublayer < tau_cone.size(); idx_sublayer++ ){
        auto & tau_sublayer = tau_cone[idx_sublayer];
        auto & sig_sublayer = sig_cone[idx_sublayer];
        std::vector<size_t> pos_isect;
        std::set_intersection(tau_sublayer.begin(),tau_sublayer.end(),
                              sig_sublayer.begin(),sig_sublayer.end(),
                              back_inserter(pos_isect));
        res_cone.emplace_back(pos_isect);
    }

    for(const auto & [i,c] : iter::enumerate_reverse(tau_cone)){
        fmt::print("tau{}: {}\n",i,c );
    }
    for(const auto & [i,c] : iter::enumerate_reverse(sig_cone)){
        fmt::print("sig{}: {}\n",i,c );
    }
    for(const auto & [i,c] : iter::enumerate_reverse(res_cone)){
        fmt::print("res{}: {}\n",i,c );
    }
    return res_cone;
}


qm::Scalar qm::lbit::get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> & unitary_layers, const Eigen::Matrix2cd & tau, size_t pos_tau, const Eigen::Matrix2cd & sig, size_t pos_sig){


    // Generate gates for the operators
    auto tau_gate = qm::Gate{tau, {pos_tau}, {2}};
    auto sig_gate = qm::Gate{sig, {pos_sig}, {2}};

    auto g = tau_gate;
    tools::log->info("Computing Trace (tau_{} sig_{})", pos_tau, pos_sig);
    auto lc2 = get_pos_in_lightcone2(unitary_layers,pos_tau, pos_sig);

    std::vector<std::string> net;
    std::vector<std::string> log;
    std::string empty_layer;
    size_t uw = 7; // Width of a unitary box
    size_t hw = 3; // Half-width of a unitary box
    size_t op = 3; // Overlap
    for(const auto & [idx_layer,layer] : iter::enumerate(unitary_layers)){
        // Generate
        if(g.pos.empty()) break;
        if(layer.empty()) continue;
        std::vector<size_t> gate_sequence;
        size_t gate_size = layer.front().pos.size();
        size_t pos_max   = layer.back().pos.back();
        for(size_t offset = 0; offset < gate_size; offset++) {
            if(offset + gate_size > pos_max +1) break;
            auto off_idx = num::range<size_t>(offset, pos_max - gate_size + 2, gate_size);
            if(num::mod<size_t>(offset, 2) == 1) std::reverse(off_idx.begin(), off_idx.end()); // If odd, reverse the sequence
            gate_sequence.insert(gate_sequence.end(), off_idx.begin(), off_idx.end());
        }
        tools::log->info("Gate sequence: {}", gate_sequence);
        if(net.empty()){
            empty_layer = fmt::format("{0:^{1}}"," ", pos_max*(uw-op) + op);
            net.emplace_back(empty_layer);
            net.back().replace(pos_tau * (uw-op), hw, fmt::format("[{1:^{0}}]",hw-2, fmt::format("{}",pos_tau)));
            log.emplace_back(fmt::format("insert tau [{}] -> now {}", pos_tau, g.pos));
        }
        std::vector<std::string> layer_str(gate_size, empty_layer);
        std::vector<std::string> story_str(gate_size);
        for(const auto &[idx, pos_gate] : iter::enumerate(gate_sequence)) {
            if(g.pos.empty()) break;
            auto &u = layer.at(pos_gate);
            auto idx_sublayer = num::mod<size_t>(pos_gate, gate_size);
            // Going through the sequence first forward, then backward, we are handed u gates which may or may not connect to our current g gate.
            // For a successful connection, at least one pos in g should be present in u gate.
            // After connecting, trace away legs of g that are outside of the light-cone intersection between tau and sigma.
            // Always trace the position furthest away from the target operator


            // Check if g.pos and u.pos have sites in common
            std::vector<size_t> pos_isect;
            std::set_intersection(g.pos.begin(),g.pos.end(),
                                  u.pos.begin(),u.pos.end(),
                                  back_inserter(pos_isect));

            if(not pos_isect.empty()){
                // Found a matching u. Connect it
                auto pos_old = g.pos;
                g = g.insert(u);
                tools::log->info("insert: layer [{},{}] = {} | u.pos {} | g.pos {} -> {} | intersect {}",idx_layer,idx_sublayer,2*idx_layer+idx_sublayer, u.pos, pos_old, g.pos,pos_isect);
                layer_str[idx_sublayer].replace(u.pos.front() * (uw-op), uw, fmt::format("[{1:^{0}}]",uw-2, fmt::format("{:<2},{:>2}",u.pos.front(), u.pos.back())));
                story_str[idx_sublayer].append(fmt::format("insert u{} -> now {} ", u.pos, g.pos));
            }

            // Check if g.pos has sites outside of the light-cone intersection
//            std::vector<size_t> pos_allowed = get_pos_in_lightcone(gate_size, pos_max, unitary_layers.size(), idx_layer, idx_sublayer, pos_tau, pos_sig);
            std::vector<size_t> pos_allowed = lc2[2*idx_layer+idx_sublayer];
            std::vector<size_t> pos_outside;
            std::set_difference(g.pos.begin(),g.pos.end(),
                                pos_allowed.begin(), pos_allowed.end(),
                                back_inserter(pos_outside));
            if(not pos_outside.empty()){
                // Found positions outside of the light cone. Trace them
                auto pos_old = g.pos;
                g = g.trace_pos(pos_outside);
                tools::log->info("trace : layer [{},{}] = {} | u.pos {} | g.pos {} -> {} | allowed {} | outside {}",idx_layer,idx_sublayer,2*idx_layer+idx_sublayer,u.pos, pos_old, g.pos, pos_allowed, pos_outside);
                story_str[idx_sublayer].append(fmt::format("trace {} -> now {} ", pos_outside, g.pos));
            }


        }
        for(auto & l : layer_str) net.emplace_back(l);
        for(auto & s : story_str) log.emplace_back(s);
    }

    Scalar result;

    if(g.pos.empty()){
        long distance = std::abs(static_cast<long>(pos_tau) - static_cast<long>(pos_sig));
        long max_dist = static_cast<long>(unitary_layers.size()*unitary_layers.front().front().pos.size());
        if(distance <= max_dist)
            tools::log->warn("Expected last gate to have size >= 1, since |pos_tau {} - pos_sig {}| = {} <= {}. Got size {}", pos_tau, pos_sig, distance, max_dist,g.pos.size()  );
//            throw std::runtime_error(fmt::format("Expected last gate to have size >= 1, since |pos_tau {} - pos_sig {}| = {} <= {}. Got size {}", pos_tau, pos_sig, distance, max_dist,g.pos.size() ));
        if(g.op.dimension(0) * g.op.dimension(1) != 1) throw std::runtime_error(fmt::format("Expected empty gate to have scalar op: Got dims {}", g.op.dimensions()));
        log.emplace_back(fmt::format("-> result = {:.8f}{:+.8f}i", result.real(),result.imag()));
        result = g.op.coeff(0);
    }else{
        //    if(g.pos.empty()) throw std::logic_error(fmt::format("Expected last gate to have size >= 1. Got size {}", g.pos.size()));
        // In the last step we connect the sigma operator and trace everything down to a scalar
        result = g.connect_under(sig_gate).trace();
        net.emplace_back(empty_layer);
        net.back().replace(pos_sig * (uw-op), hw, fmt::format("[{1:^{0}}]",hw-2, fmt::format("{}",pos_sig)));
        log.emplace_back(fmt::format("insert sigma [{}] -> now {} -> result = {:.8f}{:+.8f}i", pos_sig, g.pos,result.real(),result.imag()));
    }


    tools::log->info("Computed Trace (tau_{} sig_{})", pos_tau, pos_sig);
    for(const auto & [idx,layer] :iter::enumerate_reverse(net)){
        if(idx == 0)
            std::cout << fmt::format("tau  :") << layer << " | log: " << log[idx] << std::endl;
        else if( idx == net.size()-1)
            std::cout << fmt::format("sig  :") << layer << " | log: " << log[idx] << std::endl;
        else
            std::cout << fmt::format("u[{:^2}]:",idx-1) << layer << " | log: " << log[idx] << std::endl;
    }
    return result;

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
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO(1, 1, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);

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
    auto                     I        = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO(2, 2, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);

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
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_I(1, 1, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_S(1, 1, spin_dim, spin_dim);
    MPO_I.setZero();
    MPO_S.setZero();
    MPO_I.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);
    MPO_S.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(Eigen::MatrixXcd::Identity(spin_dim, spin_dim));

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
    auto                     I        = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_S(2, 2, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_I(2, 2, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_P(2, 2, spin_dim, spin_dim);
    MPO_S.setZero();
    MPO_I.setZero();
    MPO_P.setZero();

    MPO_S.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix1);
    MPO_S.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix2);
    MPO_I.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO_I.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO_P.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO_P.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix1);

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
    auto                     I          = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_I(num_paulis, num_paulis, spin_dim, spin_dim);
    MPO_S.setZero();
    MPO_I.setZero();
    for(long diag_pos = 0; diag_pos < num_paulis; diag_pos++) {
        MPO_S.slice(Eigen::array<long, 4>{diag_pos, diag_pos, 0, 0}, extent4).reshape(extent2) =
            Textra::MatrixTensorMap(paulimatrices[static_cast<size_t>(diag_pos)]);
        MPO_I.slice(Eigen::array<long, 4>{diag_pos, diag_pos, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
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
    long                  spin_dim = 2;
    Eigen::array<long, 4> extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::array<long, 4> offset4  = {0, 0, 0, 0};

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
                mpo.slice(offset4, extent4).reshape(extent2) = Textra::MatrixToTensor(pauli);
                break;
            }
            case RandomizerMode::SHUFFLE: {
                std::shuffle(pauli_idx.begin(), pauli_idx.end(), rnd::internal::rng);
                [[fallthrough]];
            }
            case RandomizerMode::ASIS: {
                auto num_paulis = static_cast<long>(paulimatrices.size());
                mpo.resize(num_paulis, num_paulis, spin_dim, spin_dim);
                for(long idx = 0; idx < num_paulis; idx++) {
                    auto        uidx                             = static_cast<size_t>(idx);
                    const auto &pauli                            = paulimatrices[pauli_idx[uidx]];
                    offset4                                      = {idx, idx, 0, 0};
                    mpo.slice(offset4, extent4).reshape(extent2) = Textra::MatrixToTensor(pauli);
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
    long                  num_paulis = static_cast<long>(paulimatrices.size());
    long                  spin_dim   = 2;
    Eigen::array<long, 4> extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */

    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    for(size_t site = 0; site < sites; site++) {
        Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
        MPO_S.setZero();
        for(long idx = 0; idx < num_paulis; idx++) {
            auto        uidx                               = static_cast<size_t>(idx);
            auto        coeff                              = 1 + rnd::uniform_double_box(uniform_dist_widths[uidx]);
            auto        offset4                            = Eigen::array<long, 4>{idx, idx, 0, 0};
            const auto &pauli                              = paulimatrices[uidx];
            MPO_S.slice(offset4, extent4).reshape(extent2) = Textra::MatrixToTensor(coeff * pauli);
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
