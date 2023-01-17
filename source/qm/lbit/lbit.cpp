#include "../lbit.h"
#include "../spin.h"
#include "config/debug.h"
#include "config/enums.h"
// #include "config/settings.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "io/fmt.h"
#include "io/spdlog.h"
#include "math/fit.h"
#include "math/linalg.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/stat.h"
#include "math/tenx.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include <algorithm>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

namespace settings {
    inline constexpr bool debug_circuit = false;
}

using cplx = qm::cplx;

std::vector<qm::Gate> qm::lbit::get_unitary_2gate_layer(size_t sites, double fmix) {
    /*! Returns a set of unitary two site operators used to transform between physical and l-bit representations
     *
     *
     @verbatim
                   0      1                0
                   |      |                |
                 [ exp(-ifM) ]  ---> [ exp(-ifM) ]
                   |      |               |
                   2      3               1

     @endverbatim

     Where M is defined as
        θ³n[i]n[i+1] +
        θ²n[i](1-n[i+1]) +
        θ¹(1-n[i])n[i+1] +
        θ⁰(1-n[i])(1-n[i+1]) +
        c  σ+σ- +
        c* σ-σ+

     and n[i] = 0.5(σz[i] + 1) acts on site i.
     Alternatively, we could write this in terms of σz operators, as
        θ³/4 (1+σz[i])(1+σz[i+1]) +
        θ²/4 (1+σz[i])(1-σz[i+1]) +
        θ¹/4 (1-σz[i])(1+σz[i+1]) +
        θ⁰/4 (1-σz[i])(1-σz[i+1]) +
        c  σ+σ- +
        c* σ-σ+

    */

    tools::log->trace("Generating twosite unitaries");
    auto                  SZ        = qm::spin::half::gen_twobody_spins(qm::spin::half::sz); // We use these as matrices
    auto                  SP        = qm::spin::half::gen_twobody_spins(qm::spin::half::sp); // We use these as matrices
    auto                  SM        = qm::spin::half::gen_twobody_spins(qm::spin::half::sm); // We use these as matrices
    auto                  ID        = qm::spin::half::gen_twobody_spins(qm::spin::half::id); // We use these as matrices
    auto                  N         = std::vector<Eigen::Matrix4cd>{0.5 * (ID[0] + SZ[0]), 0.5 * (ID[1] + SZ[1])};
    auto                  spin_dims = std::vector<long>{2l, 2l};
    std::vector<qm::Gate> unitaries;
    unitaries.reserve(sites - 1);
    for(size_t idx = 0; idx < sites - 1; idx++) {
        double               th0 = rnd::normal(0, 1);
        double               th1 = rnd::normal(0, 1);
        double               th2 = rnd::normal(0, 1);
        double               th3 = rnd::normal(0, 1);
        std::complex<double> c(rnd::normal(0, 1), rnd::normal(0, 1));
        auto                 indices = std::vector<size_t>{idx, idx + 1};
        Eigen::Matrix4cd     H       = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
                             SP[0] * SM[1] * c + SP[1] * SM[0] * std::conj(c);

        // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
        // the kronecker product that generated two-site gates above has indexed right-to-left
        //         0                   1      0              0      1                0
        //         |                   |      |              |      |                |
        //   [ exp(-ifM) ]  --->    [ exp(-ifM) ]   --->  [ exp(-ifM) ]  --->  [ exp(-ifM) ]
        //         |                   |      |              |      |                |
        //         1                   3      2              2      3                1
        Eigen::Tensor<cplx, 2> M_shuffled = tenx::TensorMap(H, 2, 2, 2, 2).shuffle(tenx::array4{1, 0, 3, 2}).reshape(tenx::array2{4, 4});
        Eigen::MatrixXcd       expifM     = (imn * fmix * tenx::MatrixMap(M_shuffled)).exp();
        unitaries.emplace_back(tenx::TensorMap(expifM), indices, spin_dims);
    }
    if constexpr(settings::debug) {
        // Sanity check
        for(const auto &u : unitaries)
            if(not tenx::MatrixMap(u.op).isUnitary()) throw except::logic_error("u is not unitary!");
    }

    return unitaries;
}

std::vector<qm::Gate> qm::lbit::get_unitary_2gate_layer_blocked(size_t sites, double fmix, const std::vector<double> &fields,
                                                                [[maybe_unused]] double fieldvar) {
    /*! Returns a set of unitary two site operators used to transform between physical and l-bit representations
     *
     *
     @verbatim
                   0      1                0
                   |      |                |
                 [ exp(-ifM) ]  ---> [ exp(-ifM) ]
                   |      |               |
                   2      3               1
     @endverbatim

     Where M is defined as
        θ³n[i]n[i+1] +
        θ²n[i](1-n[i+1]) +
        θ¹(1-n[i])n[i+1] +
        θ⁰(1-n[i])(1-n[i+1]) +
        c  σ+σ- +
        c* σ-σ+

     and n[i] = 0.5(σz[i] + 1) acts on site i.
     Alternatively, we could write this in terms of σz operators, as
        θ³/4 (1+σz[i])(1+σz[i+1]) +
        θ²/4 (1+σz[i])(1-σz[i+1]) +
        θ¹/4 (1-σz[i])(1+σz[i+1]) +
        θ⁰/4 (1-σz[i])(1-σz[i+1]) +
        c  σ+σ- +
        c* σ-σ+
    */

    tools::log->trace("Generating correlated twosite unitaries");
    if(fields.size() != sites) throw except::logic_error("fields.size() {} != sites {}", fields.size(), sites);
    constexpr bool        kroneckerSwap = false;
    auto                  SZ            = qm::spin::half::gen_twobody_spins(qm::spin::half::sz, kroneckerSwap); // We use these as matrices
    auto                  SP            = qm::spin::half::gen_twobody_spins(qm::spin::half::sp, kroneckerSwap); // We use these as matrices
    auto                  SM            = qm::spin::half::gen_twobody_spins(qm::spin::half::sm, kroneckerSwap); // We use these as matrices
    auto                  ID            = qm::spin::half::gen_twobody_spins(qm::spin::half::id, kroneckerSwap); // We use these as matrices
    auto                  N             = std::vector<Eigen::Matrix4cd>{0.5 * (ID[0] + SZ[0]), 0.5 * (ID[1] + SZ[1])};
    auto                  spin_dims     = std::vector<long>{2l, 2l};
    std::vector<qm::Gate> unitaries;
    unitaries.reserve(sites - 1);
    for(size_t idx = 0; idx < sites - 1; idx++) {
        // This 2-site gate connects sites idx and idx+1
        double diff = std::abs(fields[idx] - fields[idx + 1]);
        double dcay = std::exp(-2.0 * diff); // Squared exponential decay
        auto   th0  = rnd::normal(0, 1);
        auto   th1  = rnd::normal(0, 1);
        auto   th2  = rnd::normal(0, 1);
        auto   th3  = rnd::normal(0, 1);
#pragma message "Trying c with larger variance = 2²"
        auto c = std::complex<double>(rnd::normal(0, 2), rnd::normal(0, 2)) * dcay; // complex normal random variable

        auto             indices = std::vector<size_t>{idx, idx + 1};
        Eigen::Matrix4cd M       = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
                             SP[0] * SM[1] * c + SP[1] * SM[0] * std::conj(c);

        // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
        // the kronecker product that generated two-site gates above has indexed right-to-left
        //         0                   1      0              0      1                0
        //         |                   |      |              |      |                |
        //   [ exp(-ifM) ]  --->    [ exp(-ifM) ]   --->  [ exp(-ifM) ]  --->  [ exp(-ifM) ]
        //         |                   |      |              |      |                |
        //         1                   3      2              2      3                1
        Eigen::Tensor<cplx, 2> M_shuffled = tenx::TensorMap(M, 2, 2, 2, 2).shuffle(tenx::array4{1, 0, 3, 2}).reshape(tenx::array2{4, 4});
        Eigen::MatrixXcd       expifM     = (imn * fmix * tenx::MatrixMap(M_shuffled)).exp();
        unitaries.emplace_back(tenx::TensorMap(expifM), indices, spin_dims);
    }
    if constexpr(settings::debug) {
        // Sanity check
        for(const auto &u : unitaries)
            if(not tenx::MatrixMap(u.op).isUnitary()) throw except::logic_error("u is not unitary!");
    }

    return unitaries;
}

Eigen::Tensor<cplx, 2> qm::lbit::get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<cplx, 2> &H) {
    // Given a matrix H, this returns exp(delta_t * H)
    // For time evolution, just make sure delta_t = -i*d,  where d is a (small) real positive number.
    return tenx::TensorCast((delta_t * tenx::MatrixMap(H)).exp());
}

std::vector<qm::Gate> qm::lbit::get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite, double id_threshold) {
    std::vector<Gate> time_evolution_gates;
    time_evolution_gates.reserve(hams_nsite.size());
    for(auto &h : hams_nsite) {
        time_evolution_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
        if(tenx::isIdentity(time_evolution_gates.back().op, id_threshold)) {
            tools::log->trace("get_time_evolution_gates: ignoring time evolution swap gate {} == I +- {:.2e}", time_evolution_gates.back().pos, id_threshold);
            time_evolution_gates.pop_back(); // Skip this gate if it is just an identity.
        }
    }
    if constexpr(settings::debug) {
        for(auto &t : time_evolution_gates)
            if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                throw except::runtime_error("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op));
            }
    }
    time_evolution_gates.shrink_to_fit();
    return time_evolution_gates;
}

std::vector<qm::SwapGate> qm::lbit::get_time_evolution_swap_gates(cplx delta_t, const std::vector<qm::SwapGate> &hams_nsite, double id_threshold) {
    std::vector<SwapGate> time_evolution_swap_gates;
    time_evolution_swap_gates.reserve(hams_nsite.size());
    size_t count_ignored = 0;
    for(auto &h : hams_nsite) {
        time_evolution_swap_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
        if(tenx::isIdentity(time_evolution_swap_gates.back().op, id_threshold)) {
            count_ignored++;
            tools::log->trace("get_time_evolution_swap_gates: ignoring time evolution swap gate {} == I +- {:.2e}", time_evolution_swap_gates.back().pos,
                              id_threshold);
            time_evolution_swap_gates.pop_back(); // Skip this gate if it is just an identity.
        } else {
            time_evolution_swap_gates.back().generate_swap_sequences();
        }
    }
    if constexpr(settings::debug) {
        for(auto &t : time_evolution_swap_gates)
            if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                throw except::runtime_error("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op));
            }
    }
    time_evolution_swap_gates.shrink_to_fit();
    tools::log->debug("get_time_evolution_swap_gates: ignored {} time evolution swap gates: I +- {:.2e}", count_ignored, id_threshold);
    return time_evolution_swap_gates;
}

std::vector<Eigen::Tensor<cplx, 2>> qm::lbit::get_time_evolution_operators_2site(size_t sites, cplx delta_t,
                                                                                 const std::vector<Eigen::Tensor<cplx, 2>> &hams_2site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{j,j+1} and h_{j+1,j+2} commute. Therefore we can immediately use the relation
    //      exp(-i*dt *[h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}]) =  exp(-i*dt [h_{j,j+1}]) * exp(-i*dt*[h_{j+1,j+2}]) * ... * exp(-i*dt*[h_{L-2, L-1}])
    // without passing through the Suzuki-Trotter decomposition.
    // Here we expect "hams_2site" to contain terms like  h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}.

    if(hams_2site.size() != sites - 1) throw except::logic_error("Wrong number of twosite hamiltonians: {}. Expected {}", hams_2site.size(), sites - 1);

    std::vector<Eigen::Tensor<cplx, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_2site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

std::vector<Eigen::Tensor<cplx, 2>> qm::lbit::get_time_evolution_operators_3site(size_t sites, cplx delta_t,
                                                                                 const std::vector<Eigen::Tensor<cplx, 2>> &hams_3site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{i,j,k} and h_{l,m,n} commute. Therefore we can immediately use the relation
    // exp(A + B) = exp(A)exp(B)
    // without passing through the Suzuki-Trotter decomposition.

    if(hams_3site.size() != sites - 2) throw except::logic_error("Wrong number of three-site hamiltonians: {}. Expected {}", hams_3site.size(), sites - 2);

    std::vector<Eigen::Tensor<cplx, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_3site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

qm::cplx qm::lbit::get_lbit_exp_value3(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi,
                                       const Eigen::Matrix2cd &szj, size_t pos_szj, long sites) {
    /*! \brief Calculates the operator overlap O(i,j) = Tr(𝜌'_i σ^z_j) / Tr(𝜌'_i) = Tr(τ^z_i σ^z_j) / 2^L

        Consider first the density matrix

             ρ_i = 1/2^L  (σ^z_i + 1 ) ⊗ I_i',

        where I_i' is the identity matrix on sites j != i.
        Then ρ_i represents a state |psi><psi| where site i is at level 0 with probability 1,
        and the others sites are in a cat-state of level 0 or 1. E.g. ρ could be a state with
        magnetization 1 at site i, and 0 elsewhere, or ρ could be a state with 1 particle at
        site i, and cat state of 0 and 1 particles elsewhere.

        Applying the unitary transformation does not change the trace, since it is just a
        change of basis. We get

             ρ'_i = U† ρ_i U
                  = 1/2^L (U† σ^z_i U + 1)
                  = 1/2^L (τ^z_i + 1) ,

        where the l-bit τ^z_i = U† σ^z_i U acts non-trivially on all sites. Carrying out the
        trace gives us

            Tr(ρ'_i σ^z_j) / Tr(ρ'_i) = 1/2^L Tr(τ^z_i σ^z_j).

        where we used Tr(1) = 2^L, Tr(ρ_i) = Tr(ρ'_i) = 1 and Tr(τ^z_i) = 0.

        We expect an l-bit fully localized at site i to give O(i,i) = 1,  and a fully
        delocalized l-bit to give O(i,j) = 0.5.

        Note that when the light-cone from site i can't reach site j in the unitary circuit, then
        τ^z_i has no support where σ^z_j connects.  Therefore, we get

               Tr(τ^z_i σ^z_j) = Tr(τ^z_i ⊗ σ^z_j) = Tr(τ^z_i) Tr(σ^z_j) = 0

        Read more here:
        https://link.aps.org/doi/10.1103/PhysRevB.91.085425
        https://onlinelibrary.wiley.com/doi/10.1002/andp.201600322
    */

    // Generate gates for the operators
    if constexpr(settings::debug_circuit) tools::log->trace("Computing Tr (τ_{} σ_{}) / Tr(σ_{})", pos_szi, pos_szj, pos_szj);
    auto t_olap = tid::ur();
    t_olap.tic();

    auto result          = cplx(0, 0);
    auto szi_gate        = qm::Gate(szi, {pos_szi}, {2l}); //
    auto szj_gate        = qm::Gate(szj, {pos_szj}, {2l});
    auto intersection    = qm::get_lightcone_intersection(unitary_layers, pos_szi, pos_szj);
    auto unitary_slayers = qm::get_lightcone_gate_selection(unitary_layers, intersection, false); // Selected gates in each layer
    auto is_disconnected = std::any_of(intersection.begin(), intersection.end(), [](auto &layer) { return layer.empty(); });
    if(is_disconnected) {
        if constexpr(settings::debug_circuit) tools::log->trace("σzi:{} and σzj:{} are disconnected -> result = {:.6f}", pos_szi, pos_szj, result);
        return result;
    }
    // Setup debug printing of unitary circuits
    size_t                  uw        = 7;                                // Width of a unitary 2-site gate box
    size_t                  op        = 3;                                // Overlap of a unitary 2-site gate box
    size_t                  tw        = 6;                                // Tag width
    size_t                  mp        = static_cast<size_t>(sites) - 1ul; // Max pos
    auto                    empty_str = fmt::format("{0:^{1}}", " ", tw + mp * (uw - op) + op);
    std::deque<std::string> net(2 + unitary_slayers.size(), empty_str); // Great for debugging
    std::deque<std::string> log(2 + unitary_slayers.size());            // Great for debugging

    auto all_empty = [](const auto &slayers) {
        return std::all_of(slayers.begin(), slayers.end(), [](const auto &slayer) -> bool { return slayer.empty(); });
    };
    auto num_gates_L = [](const auto &slayers, size_t idx = 0) {
        size_t num = 0;
        for(size_t i = idx + 1; i < slayers.size(); i++) {
            if(slayers.at(i - 1).empty()) return num;
            if(slayers.at(i).empty()) return num;
            size_t posLL = slayers.at(i).front().pos.front();
            size_t posL  = slayers.at(i - 1).front().pos.front();
            if(posLL + 1 == posL)
                num++;
            else
                break;
        }
        return num;
    };
    auto num_gates_R = [](const auto &slayers, size_t idx = 0) {
        size_t num = 0;
        for(size_t i = idx + 1; i < slayers.size(); i++) {
            if(slayers.at(i - 1).empty()) return num;
            if(slayers.at(i).empty()) return num;
            size_t posRR = slayers.at(i).back().pos.back();
            size_t posR  = slayers.at(i - 1).back().pos.back();
            if(posR + 1 == posRR)
                num++;
            else
                break;
        }
        return num;
    };
    auto slayer_width = [&]() -> size_t {
        size_t width = 0;
        for(const auto &s : unitary_slayers) {
            size_t w = 0;
            for(const auto &l : s) w += l.pos.size();
            width = std::max(width, w);
        }
        return width;
    };
    auto slayer_diagl = [&]() -> size_t {
        size_t lb = -1ul; // left bottom position
        size_t rb = -1ul; // right bottom position
        size_t lt = -1ul; // left  top position
        size_t rt = -1ul; // right top position
        for(const auto &us : unitary_slayers) {
            if(not us.empty()) {
                if(lb == -1ul) lb = us.front().pos.front();
                if(rb == -1ul) rb = us.front().pos.back();
                lb = std::min(lb, us.front().pos.front());
            }
        }
        for(const auto &us : iter::reverse(unitary_slayers)) {
            if(not us.empty()) {
                if(rt == -1ul) rt = us.front().pos.back();
                if(lt == -1ul) lt = us.front().pos.front();
                rt = std::max(rt, us.back().pos.back());
            }
        }
        return std::max(rb - lb, rt - lt) + 1;
    };

    auto slayer_diagr = [&]() -> size_t {
        size_t lb = -1ul; // left bottom position
        size_t rb = -1ul; // right bottom position
        size_t lt = -1ul; // left  top position
        size_t rt = -1ul; // right top position
        for(const auto &us : unitary_slayers) {
            if(not us.empty()) {
                if(lb == -1ul) lb = us.back().pos.front();
                if(rb == -1ul) rb = us.back().pos.back();
                rb = std::max(rb, us.back().pos.back());
            }
        }
        for(const auto &us : iter::reverse(unitary_slayers)) {
            if(not us.empty()) {
                if(rt == -1ul) rt = us.back().pos.back();
                if(lt == -1ul) lt = us.back().pos.front();
                lt = std::min(lt, us.front().pos.front());
            }
        }
        return std::max(rb - lb, rt - lt) + 1;
    };

    auto width = slayer_width();
    auto diagr = slayer_diagr();
    auto diagl = slayer_diagl();

    bool go_diagonal = std::min(diagl, diagr) < width; // Avoid long thin diagonals.

    // Start contracting selected unitary gates bottom to top.
    auto g = szi_gate; // The gate that accumulates everything. Starts at σzi position
    if constexpr(settings::debug_circuit) {
        szi_gate.draw_pos(net.front(), "szi  :");
        log.front().append(fmt::format("insert σzi{}", szi_gate.pos));
    }
    if(go_diagonal) {
        // Decide to take the left or right diagonal (whatever is cheapest)
        bool go_right = diagr <= diagl;
        while(not all_empty(unitary_slayers)) {
            for(auto &&[idx_slayer, slayer] : iter::enumerate(unitary_slayers)) {
                if(slayer.empty()) continue; // Already applied the whole slayer
                auto &layer_str = net.at(idx_slayer + 1);
                auto &story_str = log.at(idx_slayer + 1);
                auto  numL      = num_gates_L(unitary_slayers, idx_slayer); // number of gates remaining to take going left
                auto  numR      = num_gates_R(unitary_slayers, idx_slayer); // number of gates remaining to take going right

                std::vector<size_t> pos_out;
                if(go_right) {
                    //                    tools::log->trace("-> insert u[{}]:{}", idx_slayer, slayer.back().pos);
                    if constexpr(settings::debug_circuit) {
                        slayer.back().draw_pos(layer_str, fmt::format("u[{:^2}]:", idx_slayer));
                        story_str.append(fmt::format("insert u{} ", slayer.back().pos));
                    }
                    g       = g.insert(slayer.back());
                    pos_out = slayer.back().pos_difference(intersection.at(idx_slayer + 1));
                    slayer.pop_back();
                } else {
                    //                    tools::log->trace("<- insert u[{}]:{}", idx_slayer, slayer.front().pos);
                    if constexpr(settings::debug_circuit) {
                        slayer.front().draw_pos(layer_str, fmt::format("u[{:^2}]:", idx_slayer));
                        story_str.append(fmt::format("insert u{} ", slayer.front().pos));
                    }
                    g       = g.insert(slayer.front());
                    pos_out = slayer.front().pos_difference(intersection.at(idx_slayer + 1));
                    slayer.pop_front();
                }

                // Collect positions that could be traced
                if(not pos_out.empty()) {                                                      // Trace positions outside the light cone intersection
                    pos_out.erase(std::unique(pos_out.begin(), pos_out.end()), pos_out.end()); // Keep unique elements
                    //                    tools::log->trace("trace[{}]:{}", idx_slayer, pos_out);
                    g    = g.trace_pos(pos_out);
                    g.op = g.op / g.op.constant(std::pow(2, pos_out.size())); // Normalize by dividing the trace of each 2x2 identity.
                    if constexpr(settings::debug_circuit) story_str.append(fmt::format("trace{} ", pos_out));
                }

                if((numR == 0 and go_right) or (numL == 0 and not go_right)) {
                    if constexpr(settings::debug_circuit) story_str.append(fmt::format("break {} ", g.pos));
                    break; // Start from the bottom of the circuit
                }
            }
        }
    } else {
        for(auto &&[idx_slayer, slayer] : iter::enumerate(unitary_slayers)) {
            auto &layer_str = net.at(idx_slayer + 1);
            auto &story_str = log.at(idx_slayer + 1);
            for(auto &sgate : slayer) {
                if(g.has_pos(sgate.pos)) {
                    g = g.insert(sgate);
                    //                    tools::log->trace("insert u[{}]:{}", idx_slayer, sgate.pos);
                    if constexpr(settings::debug_circuit) {
                        sgate.draw_pos(layer_str, fmt::format("u[{:^2}]:", idx_slayer));
                        story_str.append(fmt::format("insert u{} ", sgate.pos));
                    }
                    auto pos_out = sgate.pos_difference(intersection.at(idx_slayer + 1));
                    if(not pos_out.empty()) {
                        // Trace positions outside the light cone intersection
                        g    = g.trace_pos(pos_out);
                        g.op = g.op / g.op.constant(std::pow(2, pos_out.size())); // Normalize by dividing the trace of each 2x2 identity.
                        if constexpr(settings::debug_circuit) story_str.append(fmt::format("trace{} ", pos_out));
                    }
                }
            }
            if constexpr(settings::debug_circuit) story_str.append(fmt::format("now{} ", g.pos));
        }
    }

    if(g.has_pos(szj_gate.pos)) {
        // Connect σ^z_j at the top
        szj_gate.draw_pos(net.back(), "szj  :");
        log.back().append(fmt::format("insert σzj{} ", szj_gate.pos));
        g = szj_gate.connect_above(g);
    }
    // In the last step we trace everything down to a cplx
    result = g.trace();
    result /= std::pow(2, g.pos.size()); // Normalize by dividing the trace of each 2x2 identity.
    t_olap.toc();
    if constexpr(settings::debug_circuit) {
        log.back().append(fmt::format("result = {:.6f} | time {:.3e}", result, t_olap.get_time()));
        for(const auto &[idx, layer] : iter::enumerate_reverse(net)) tools::log->debug("{} | {}", layer, log[idx]);
    }
    return result;
}

Eigen::Tensor<qm::cplx, 2> qm::lbit::get_lbit_support(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites) {
    auto                       ssites = static_cast<long>(sites);
    Eigen::Tensor<qm::cplx, 2> lbit_overlap(ssites, ssites);

    /*! \brief Calculates the operator overlap O(i,j) = Tr(𝜌_i σ^z_j) / Tr(𝜌_j)
                                                      = Tr(τ^z_i σ^z_j) / 2^L
                                                      = Tr(U†σ^z_iU σ^z_j) / 2^L
        Where
            𝜌_i   = 1/2^L (1 + τ^z_i)   a density matrix describing a state localized around site i
            τ^z_i = U† σ^z_i  U,        is an l-bit operator at site i
            U                           is a unitary transformation in the form of a finite-depth circuit.
        The second equality comes from the fact that Tr(τ^z_i) = 0 (l-bit operators are traceless) and Tr(1) = 2^L.
        An l-bit fully localized at site i gives O(i,i) = 1,  and when fully delocalized one gets O(i,j) = 0.5

        Read more here:
        https://link.aps.org/doi/10.1103/PhysRevB.91.085425
        https://onlinelibrary.wiley.com/doi/10.1002/andp.201600322
    */
#pragma omp parallel for collapse(2) schedule(guided, 4)
    for(long j = 0; j < ssites; j++) {
        for(long i = 0; i < ssites; i++) {
            lbit_overlap(i, j) =
                qm::lbit::get_lbit_exp_value3(unitary_layers, qm::spin::half::sz, static_cast<size_t>(i), qm::spin::half::sz, static_cast<size_t>(j), ssites);
        }
    }
    // We require that lbit_overlap(i,j) has rows that sum up to 1
    auto sums_rowwise = tenx::MatrixMap(lbit_overlap).rowwise().sum();
    if(not sums_rowwise.cwiseAbs().isOnes(1e-4)) {
        tools::log->info("lbit_overlap: \n{}\nsums\n{}\n", linalg::tensor::to_string(lbit_overlap, 6), linalg::matrix::to_string(sums_rowwise, 6));
        throw except::logic_error("lbit overlap rows do not sum to one. Perhaps normalization is wrong");
    }
    return lbit_overlap;
}

Eigen::Tensor<cplx, 2> qm::lbit::get_lbit_permuted_support(const Eigen::Tensor<cplx, 2> &lbit_support) {
    // First, subtract the center position of each lbit, so we get L lbits centered around zero.
    // In practice, we make a cyclic permutation of the rows of lbit_support
    // In addition, we mirror the lbit along its vertical, so that we can average its left and right half together
    long                   rows = lbit_support.dimension(0);
    long                   cols = lbit_support.dimension(1);
    Eigen::Tensor<cplx, 2> lbit_permuted_support(rows, cols);
    lbit_permuted_support.setConstant(0);
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
                              std::real(lbit_support(i, j)), std::real(lbit_support(i, j_twin)),
                              0.5 * (std::real(lbit_support(i, j)) + std::real(lbit_support(i, j_twin))));
            if(std::abs(std::imag(lbit_support(i, j))) > 1e-12)
                tools::log->warn("lbit_overlap({},{}) has imaginary component : |Im({})| > 1e-12", i, j, lbit_support(i, j));
            lbit_permuted_support(i, j_perm) = 0.5 * (lbit_support(i, j) + lbit_support(i, j_twin));
        }
    }
    return lbit_permuted_support;
}

std::tuple<double, double, std::vector<double>, size_t> qm::lbit::get_characteristic_length_scale(const Eigen::Tensor<cplx, 2> &lbit_overlap_permuted) {
    // Average along each column to get an estimate of the lbit
    Eigen::Tensor<double, 1> lbit_overlap_perm_avg = lbit_overlap_permuted.real().mean(std::array<long, 1>{0}).abs(); // Average each column
    Eigen::Tensor<double, 1> lbit_overlap_perm_log = lbit_overlap_perm_avg.log();                                     // Make all values positive, and take log
    tools::log->info("lbit_overlap_perm_avg: \n{}\n", linalg::tensor::to_string(lbit_overlap_perm_avg.real(), 18, 20));
    tools::log->info("lbit_overlap_perm_log: \n{}\n", linalg::tensor::to_string(lbit_overlap_perm_log.real(), 18, 20));

    // Data becomes noisy if the exponential has decayed, so find a cutoff to get the slope using only the first part of the curve
    auto y      = std::vector<double>(lbit_overlap_perm_avg.data(), lbit_overlap_perm_avg.data() + lbit_overlap_perm_avg.size());
    auto v      = stat::find_last_valid_point(y);
    auto c      = std::count_if(y.begin(), y.begin() + static_cast<long>(v), [](auto &val) { return std::abs(val) > std::numeric_limits<double>::epsilon(); });
    auto x_full = num::range<double>(0, c);
    auto x_tail = num::range<double>(c / 2, c);
    auto y_full = tenx::span(lbit_overlap_perm_avg.data(), c);
    auto y_tail = tenx::span(y_full.data() + c / 2, y_full.data() + c);
    auto ylog_full                = tenx::span(lbit_overlap_perm_log.data(), c);
    auto ylog_tail                = tenx::span(ylog_full.data() + c / 2, ylog_full.data() + c);
    auto [slope_full, res_full]   = stat::slope(x_full, ylog_full);
    auto [slope_tail, res_tail]   = stat::slope(x_tail, ylog_tail);
    double cls_full               = 1.0 / std::abs(slope_full);
    double cls_tail               = 1.0 / std::abs(slope_tail);
    auto   fit_log_stretched_full = fit::log_stretched(x_full, ylog_full, {0.9, 1.0, 1.0});
    auto   fit_log_stretched_tail = fit::log_stretched(x_tail, ylog_tail, {0.9, 1.0, 1.0});
    auto   fit_log_full           = fit::log(x_full, ylog_full, {1.0, 1.0});
    auto   fit_log_tail           = fit::log(x_tail, ylog_tail, {1.0, 1.0});
    auto   fit_exp_stretched_full = fit::exp_stretched(x_full, y_full, {0.9, 1.0, 1.0});
    auto   fit_exp_stretched_tail = fit::exp_stretched(x_tail, y_tail, {0.9, 1.0, 1.0});
    auto   fit_exp_full           = fit::exp(x_full, y_full, {1.0, 1.0});
    auto   fit_exp_tail           = fit::exp(x_tail, y_tail, {1.0, 1.0});

    tools::log->info("Computed fit full slope    cls    {:>8.6f} | sse {:>8.6f} | using {} points: {::8.6f}", cls_full, res_full, ylog_full.size(), ylog_full);
    tools::log->info("Computed fit tail slope    cls    {:>8.6f} | sse {:>8.6f} | using {} points: {::8.6f}", cls_tail, res_tail, ylog_tail.size(), ylog_tail);
    tools::log->info("Computed fit full log_stretched coeffs {::>8.6f} | status {}", fit_log_stretched_full.coeffs, fit_log_stretched_full.status);
    tools::log->info("Computed fit tail log_stretched coeffs {::>8.6f} | status {}", fit_log_stretched_tail.coeffs, fit_log_stretched_tail.status);
    tools::log->info("Computed fit full log           coeffs {::>8.6f} | status {}", fit_log_full.coeffs, fit_log_full.status);
    tools::log->info("Computed fit tail log           coeffs {::>8.6f} | status {}", fit_log_tail.coeffs, fit_log_tail.status);
    tools::log->info("Computed fit full exp_stretched coeffs {::>8.6f} | status {}", fit_exp_stretched_full.coeffs, fit_exp_stretched_full.status);
    tools::log->info("Computed fit tail exp_stretched coeffs {::>8.6f} | status {}", fit_exp_stretched_tail.coeffs, fit_exp_stretched_tail.status);
    tools::log->info("Computed fit full exp           coeffs {::>8.6f} | status {}", fit_exp_full.coeffs, fit_exp_full.status);
    tools::log->info("Computed fit tail exp           coeffs {::>8.6f} | status {}", fit_exp_tail.coeffs, fit_exp_tail.status);
    return {cls_full, res_full, y, c};
}

std::tuple<double, double, std::vector<double>, size_t> qm::lbit::get_characteristic_length_scale2(const Eigen::Tensor<cplx, 2> &lbit_overlap_permuted) {
    // Average along each column to get an estimate of the lbit
    Eigen::Tensor<double, 1> lbit_overlap_perm_avg = lbit_overlap_permuted.real().abs().mean(std::array<long, 1>{0}); // Average each column
    Eigen::Tensor<double, 1> lbit_overlap_perm_log = lbit_overlap_perm_avg.log();                                     // Make all values positive, and take log
    tools::log->info("lbit_overlap_perm_avg: \n{}\n", linalg::tensor::to_string(lbit_overlap_perm_avg.real(), 18, 20));
    tools::log->info("lbit_overlap_perm_log: \n{}\n", linalg::tensor::to_string(lbit_overlap_perm_log.real(), 18, 20));

    // Data becomes noisy if the exponential has decayed, so find a cutoff to get the slope using only the first part of the curve
    auto y      = std::vector<double>(lbit_overlap_perm_avg.data(), lbit_overlap_perm_avg.data() + lbit_overlap_perm_avg.size());
    auto v      = stat::find_last_valid_point(y);
    auto c      = std::count_if(y.begin(), y.begin() + static_cast<long>(v), [](auto &val) { return std::abs(val) > std::numeric_limits<double>::epsilon(); });
    auto x_full = num::range<double>(0, c);
    auto x_tail = num::range<double>(c / 2, c);
    auto y_full = tenx::span(lbit_overlap_perm_avg.data(), c);
    auto y_tail = tenx::span(y_full.data() + c / 2, y_full.data() + c);
    auto ylog_full                = tenx::span(lbit_overlap_perm_log.data(), c);
    auto ylog_tail                = tenx::span(ylog_full.data() + c / 2, ylog_full.data() + c);
    auto [slope_full, res_full]   = stat::slope(x_full, ylog_full);
    auto [slope_tail, res_tail]   = stat::slope(x_tail, ylog_tail);
    double cls_full               = 1.0 / std::abs(slope_full);
    double cls_tail               = 1.0 / std::abs(slope_tail);
    auto   fit_log_stretched_full = fit::log_stretched(x_full, ylog_full, {0.9, 1.0, 1.0});
    auto   fit_log_stretched_tail = fit::log_stretched(x_tail, ylog_tail, {0.9, 1.0, 1.0});
    auto   fit_log_full           = fit::log(x_full, ylog_full, {1.0, 1.0});
    auto   fit_log_tail           = fit::log(x_tail, ylog_tail, {1.0, 1.0});
    auto   fit_exp_stretched_full = fit::exp_stretched(x_full, y_full, {0.9, 1.0, 1.0});
    auto   fit_exp_stretched_tail = fit::exp_stretched(x_tail, y_tail, {0.9, 1.0, 1.0});
    auto   fit_exp_full           = fit::exp(x_full, y_full, {1.0, 1.0});
    auto   fit_exp_tail           = fit::exp(x_tail, y_tail, {1.0, 1.0});

    tools::log->info("Computed fit full slope    cls    {:>8.6f} | sse {:>8.6f} | using {} points: {::8.6f}", cls_full, res_full, ylog_full.size(), ylog_full);
    tools::log->info("Computed fit tail slope    cls    {:>8.6f} | sse {:>8.6f} | using {} points: {::8.6f}", cls_tail, res_tail, ylog_tail.size(), ylog_tail);
    tools::log->info("Computed fit full log_stretched coeffs {::>8.6f} | status {}", fit_log_stretched_full.coeffs, fit_log_stretched_full.status);
    tools::log->info("Computed fit tail log_stretched coeffs {::>8.6f} | status {}", fit_log_stretched_tail.coeffs, fit_log_stretched_tail.status);
    tools::log->info("Computed fit full log           coeffs {::>8.6f} | status {}", fit_log_full.coeffs, fit_log_full.status);
    tools::log->info("Computed fit tail log           coeffs {::>8.6f} | status {}", fit_log_tail.coeffs, fit_log_tail.status);
    tools::log->info("Computed fit full exp_stretched coeffs {::>8.6f} | status {}", fit_exp_stretched_full.coeffs, fit_exp_stretched_full.status);
    tools::log->info("Computed fit tail exp_stretched coeffs {::>8.6f} | status {}", fit_exp_stretched_tail.coeffs, fit_exp_stretched_tail.status);
    tools::log->info("Computed fit full exp           coeffs {::>8.6f} | status {}", fit_exp_full.coeffs, fit_exp_full.status);
    tools::log->info("Computed fit tail exp           coeffs {::>8.6f} | status {}", fit_exp_tail.coeffs, fit_exp_tail.status);
    return {cls_full, res_full, y, c};
}

Eigen::Tensor<cplx, 2> qm::lbit::get_lbit_support_averaged(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_overlap_vec) {
    Eigen::Tensor<cplx, 2> avg;
    //    Eigen::Tensor<cplx,2> err;
    if(not lbit_overlap_vec.empty()) {
        long              rows = lbit_overlap_vec.front().dimension(0); // dim 0
        long              cols = lbit_overlap_vec.front().dimension(1); // dim 1
        size_t            reps = lbit_overlap_vec.size();               // dim 2
        std::vector<cplx> slice(reps);
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

qm::lbit::lbitSupportAnalysis qm::lbit::get_lbit_support_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec, size_t reps,
                                                                  size_t sites, const UnitaryGateProperties &ugate_props, bool randomize_fields) {
    auto t_lbit_analysis = tid::tic_scope("lbit_analysis");

    long                rows = static_cast<long>(fmix_vec.size());
    long                cols = static_cast<long>(udepth_vec.size());
    long                repl = static_cast<long>(reps);
    long                size = static_cast<long>(sites);
    lbitSupportAnalysis lbitSA(rows, cols, repl, size);
    std::array<long, 3> offset3{}, extent3{};
    std::array<long, 5> offset5{}, extent5{};

    for(size_t uidx = 0; uidx < udepth_vec.size(); uidx++) {
        for(size_t fidx = 0; fidx < fmix_vec.size(); fidx++) {
            //    for(const auto & [uidx,udep] : iter::enumerate(udepth_vec)){

            //        for (const auto & [fidx,fmix] : iter::enumerate(fmix_vec) ){
            auto t_cls = tid::ur("lbit-cls");
            t_cls.tic();
            auto                                uldx = static_cast<long>(uidx);
            auto                                fldx = static_cast<long>(fidx);
            auto                                fmix = fmix_vec[fidx];
            auto                                udep = udepth_vec[uidx];
            std::vector<double>                 cls_vec(reps);
            std::vector<double>                 sse_vec(reps);
            std::vector<Eigen::Tensor<cplx, 2>> lbit_supp_vec(reps);
            for(size_t i = 0; i < reps; i++) {
                std::vector<std::vector<qm::Gate>> layers;
                switch(ugate_props.ugate_type) {
                    case UnitaryGateType::UNIFORM: {
                        for(size_t l = 0; l < udep; l++) layers.emplace_back(qm::lbit::get_unitary_2gate_layer(sites, fmix));
                        break;
                    }
                    case UnitaryGateType::BLOCKED: {
                        if(ugate_props.fields.empty())
                            throw except::runtime_error("Fields empty, but are required when given UnitaryGateType::FIELD_CONSTRICTION");
                        if(ugate_props.fieldvar == 0)
                            tools::log->warn("Field variance is 0. This does not make sense when given UnitaryGateType::FIELD_CONSTRICTION");
                        for(size_t l = 0; l < udep; l++) {
                            if(randomize_fields) {
                                auto ugate_props_random_fields = ugate_props;
                                for(auto &field : ugate_props_random_fields.fields) {
                                    if(settings::model::lbit::distribution == "normal")
                                        field = rnd::normal(settings::model::lbit::J1_mean, settings::model::lbit::J1_wdth);
                                    else if(settings::model::lbit::distribution == "uniform")
                                        field = rnd::uniform_double_box(settings::model::lbit::J1_mean - settings::model::lbit::J1_wdth / 2.0,
                                                                        settings::model::lbit::J1_mean + settings::model::lbit::J1_wdth / 2.0);
                                    else
                                        throw except::runtime_error("Distribution not implemented: {}", settings::model::lbit::distribution);
                                }
                                layers.emplace_back(qm::lbit::get_unitary_2gate_layer_blocked(sites, fmix, ugate_props_random_fields.fields,
                                                                                              ugate_props_random_fields.fieldvar));
                            } else {
                                layers.emplace_back(qm::lbit::get_unitary_2gate_layer_blocked(sites, fmix, ugate_props.fields, ugate_props.fieldvar));
                            }
                        }
                        break;
                    }
                }

                lbit_supp_vec[i]                     = qm::lbit::get_lbit_support(layers, sites);
                offset5                              = {static_cast<long>(fidx), static_cast<long>(uidx), static_cast<long>(i), 0, 0};
                extent5                              = {1, 1, 1, lbit_supp_vec[i].dimension(0), lbit_supp_vec[i].dimension(1)};
                lbitSA.supps.slice(offset5, extent5) = lbit_supp_vec[i].real().reshape(extent5);
                lbitSA.pupps.slice(offset5, extent5) = qm::lbit::get_lbit_permuted_support(lbit_supp_vec[i]).real().reshape(extent5);
                tools::log->info("Computed u {} | f {} | rep {} | randomize fields {} | {} | {:.3e}s", udep, fmix, i, randomize_fields,
                                 enum2sv(ugate_props.ugate_type), t_cls.restart_lap());
            }

            auto lbit_supp_avg    = qm::lbit::get_lbit_support_averaged(lbit_supp_vec);
            auto lbit_pupp_avg    = qm::lbit::get_lbit_permuted_support(lbit_supp_avg);
            auto [cls, sse, y, c] = qm::lbit::get_characteristic_length_scale(lbit_pupp_avg);
            qm::lbit::get_characteristic_length_scale2(lbit_pupp_avg);
#if defined(_OPENMP)
            tools::log->info("Computed u {} | f {:.4f} | {} | lbit cls {:>8.6f} | sse {:>8.6f} | threads {} | time {:8.3f} s | decay {:2} sites: {:8.2e}", udep,
                             fmix, enum2sv(ugate_props.ugate_type), cls, sse, omp_get_max_threads(), t_cls.get_last_interval(), c, fmt::join(y, ", "));
#else
            tools::log->info("Computed u {} | f {:.4f} | {} | lbit cls {:>8.6f} | sse {:>8.6f} | time {:8.3f} s | decay {:2} sites: {:8.2e}", udep, fmix,
                             enum2sv(ugate_props.ugate_type), cls, sse, ur_iter.get_last_interval(), c, fmt::join(y, ", "));
#endif
            lbitSA.cls_avg(fldx, uldx) = cls;
            lbitSA.sse_avg(fldx, uldx) = sse;

            offset3                              = {fldx, uldx, 0};
            extent3                              = {1, 1, static_cast<long>(y.size())};
            lbitSA.decay.slice(offset3, extent3) = Eigen::TensorMap<Eigen::Tensor<double, 3>>(y.data(), extent3);
        }
    }
    return lbitSA;

    /* Example result for L = 24 result
        u 2 | f 0.0500 | decay (5 sites) : 9.94e-01, 3.22e-03, 5.24e-06, 4.90e-09, 3.98e-12, 2.75e-18
        u 3 | f 0.0500 | decay (6 sites) : 9.91e-01, 4.79e-03, 1.13e-05, 1.61e-08, 1.87e-11, 1.13e-14, 9.24e-18
        u 4 | f 0.0500 | decay (6 sites) : 9.88e-01, 6.38e-03, 2.14e-05, 4.46e-08, 7.55e-11, 7.78e-14, 8.52e-17
        u 5 | f 0.0500 | decay (7 sites) : 9.84e-01, 8.37e-03, 3.57e-05, 9.68e-08, 1.93e-10, 2.79e-13, 3.82e-16, 3.06e-19

        u 2 | f 0.1000 | decay (5 sites) : 9.75e-01, 1.31e-02, 8.83e-05, 3.08e-07, 1.05e-09, -3.12e-18,
        u 3 | f 0.1000 | decay (7 sites) : 9.62e-01, 1.99e-02, 1.90e-04, 1.17e-06, 5.83e-09, 1.33e-11, 4.27e-14, -1.83e-18
        u 4 | f 0.1000 | decay (8 sites) : 9.52e-01, 2.48e-02, 3.22e-04, 2.67e-06, 1.74e-08, 6.49e-11, 2.82e-13, 5.23e-16, -7.32e-20
        u 5 | f 0.1000 | decay (8 sites) : 9.43e-01, 2.94e-02, 4.72e-04, 5.15e-06, 4.15e-08, 2.37e-10, 1.27e-12, 4.61e-15, 1.69e-17

        u 2 | f 0.2000 | decay (5 sites) : 9.07e-01, 4.72e-02, 1.27e-03, 1.69e-05, 2.27e-07, 1.64e-18
        u 3 | f 0.2000 | decay (7 sites) : 8.61e-01, 6.96e-02, 2.79e-03, 6.40e-05, 1.31e-06, 1.16e-08, 1.62e-10, 1.27e-18
        u 4 | f 0.2000 | decay (8 sites) : 8.16e-01, 9.10e-02, 5.16e-03, 1.78e-04, 4.98e-06, 8.27e-08, 1.44e-09, 9.02e-12, 1.19e-13, 2.57e-18
        u 5 | f 0.2000 | decay (9 sites)
       : 7.73e-01, 1.11e-01, 7.90e-03, 3.41e-04, 1.15e-05, 2.55e-07, 5.52e-09, 6.77e-11, 1.12e-12, 5.91e-15, 7.47e-17, 2.38e-18

        u 2 | f 0.3000 | decay (5 sites) : 8.06e-01, 9.51e-02, 5.95e-03, 1.85e-04, 5.82e-06, 1.39e-18
        u 3 | f 0.3000 | decay (7 sites) : 7.21e-01, 1.33e-01, 1.24e-02, 6.41e-04, 3.01e-05, 5.49e-07, 1.62e-08, 1.64e-18
        u 4 | f 0.3000 | decay (9 sites) : 6.57e-01, 1.61e-01, 1.88e-02, 1.45e-03, 8.60e-05, 3.19e-06, 1.30e-07, 2.34e-09, 6.82e-11, -1.02e-18
        u 5 | f 0.3000 | decay (11
       sites): 6.27e-01, 1.70e-01, 2.44e-02, 2.36e-03, 1.80e-04, 9.49e-06, 4.92e-07, 1.53e-08, 5.85e-10, 7.55e-12, 2.20e-13, 1.39e-19

      */
}
