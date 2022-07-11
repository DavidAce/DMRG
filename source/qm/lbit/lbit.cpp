#include "../lbit.h"
#include "../spin.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "io/fmt.h"
#include "io/spdlog.h"
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
                 [ exp(-ifH) ]  ---> [ exp(-ifH) ]
                   |      |               |
                   2      3               1
     @endverbatim
    */

    tools::log->trace("Generating twosite unitaries");
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
        double               th0 = rnd::uniform_double_box(-1, 1);
        double               th1 = rnd::uniform_double_box(-1, 1);
        double               th2 = rnd::uniform_double_box(-1, 1);
        double               th3 = rnd::uniform_double_box(-1, 1);
        std::complex<double> t(rnd::uniform_double_box(-1, 1), rnd::uniform_double_box(-1, 1));

        //        double               th0 = rnd::normal(0, 1);
        //        double               th1 = rnd::normal(0, 1);
        //        double               th2 = rnd::normal(0, 1);
        //        double               th3 = rnd::normal(0, 1);
        //        std::complex<double> t    (rnd::normal(0, 1), rnd::normal(0, 1));
        //

        auto             indices = std::vector<size_t>{idx, idx + 1};
        Eigen::Matrix4cd H       = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
                             SP[0] * SM[1] * t + SP[1] * SM[0] * std::conj(t);

        if constexpr(kroneckerSwap) {
            // Here the kronecker already has index pattern left-to-right and there is no need to shuffle

            //         0               0      1
            //         |               |      |
            //   [ exp(-ifH) ]  ==  [ exp(-ifH) ]
            //        |               |      |
            //        1               2      3

            unitaries.emplace_back(tenx::TensorMap((imn * fmix * H).exp().eval()), indices, spin_dims);
        } else {
            // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
            // the kronecker product that generated two-site gates above has indexed right-to-left
            //         0                   1      0              0      1                0
            //         |                   |      |              |      |                |
            //   [ exp(-ifH) ]  --->    [ exp(-ifH) ]   --->  [ exp(-ifH) ]  --->  [ exp(-ifH) ]
            //         |                   |      |              |      |                |
            //         1                   3      2              2      3                1
            Eigen::Tensor<cplx, 2> H_shuffled = tenx::TensorMap(H, 2, 2, 2, 2).shuffle(tenx::array4{1, 0, 3, 2}).reshape(tenx::array2{4, 4});
            Eigen::MatrixXcd       expifH     = (imn * fmix * tenx::MatrixMap(H_shuffled)).exp();
            unitaries.emplace_back(tenx::TensorMap(expifH), indices, spin_dims);
        }
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

qm::cplx qm::lbit::get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &rho, size_t pos_rho,
                                      const Eigen::Matrix2cd &sig, size_t pos_sig) {
    // We calculate the operator overlap = Tr(ρ' σ^z_j) / Tr(ρ') = 1/2^L Tr(σ' σ^z_j)
    // Define
    //      ρ = (1/2^L)  (σ^z_i + 1 ) ⊗ I_i',
    // where I_i' is the identity matrix on sites j != i
    // Then ρ represents a state |psi><psi| where site i is at level 0 with probability 1, and
    // the others sites are at level 0 or 1 with probability 1/2.
    // Eg. ρ could be a state with magnetization 1 at site i, and 0 elsewhere,
    // or ρ could be a state with 1 particle at site i, and cat state of 0 and 1 particles elsewhere.
    // Applying the unitary circuit gives
    //      ρ' = U† ρ U = (1/2^L) (U† σ^z_i U + 1),
    // where σ' = U† σ^z_i U acts non-trivially on all sites.
    // Substitution and carrying out the trace gives Tr(ρ' σ^z_j) / Tr(ρ') = 1/2^L Tr(σ' σ^z_j).
    //
    // Note that when the light-cone from site i can't reach site j in the unitary circuit, then
    // ρ'_i has no support where σ^z_j connects.  Therefore we effectively get
    //        Tr(ρ'_i σ^z_j) =  Tr(ρ_i ⊗ σ^z_j) =  Tr(ρ_i) Tr(σ^z_j) = 0
    //   or alternatively
    //        Tr(ρ'_i σ^z_j) = (1/2^L) Tr([(σ^z_i +1) ⊗ I_i' ] σ^z_j) = (1/2^L) Tr( [2^L] σ^z_i ⊗ σ^z_j ) = Tr(σ^z_i) Tr(σ^z_j) = 0

    // See more about this here: https://link.aps.org/doi/10.1103/PhysRevB.91.085425

    // Generate gates for the operators
    tools::log->trace("Computing Tr (ρ_{} σ_{}) / Tr(ρ_{})", pos_rho, pos_sig);
    //    auto t_lexp   = tid::tic_scope("lbit_exp_value");
    auto rho_gate = qm::Gate(rho, {pos_rho}, {2l});
    auto sig_gate = qm::Gate(sig, {pos_sig}, {2l});
    auto g        = rho_gate; // Start with the bottom rho gate
    auto lc2      = qm::get_lightcone_intersection(unitary_layers, pos_rho, pos_sig);

    // Setup debug printing of unitary circuits
    bool                     deb = tools::log->level() <= spdlog::level::debug;
    std::vector<std::string> net; // Great for debugging
    std::vector<std::string> log; // Great for debugging
    std::string              empty_layer;
    size_t                   uw = 7; // Width of a unitary 2-site gate box
    size_t                   op = 3; // Overlap of a unitary 2-site gate box
    size_t                   tw = 6; // Tag width

    // Start contracting the unitary circuit
    for(const auto &[idx_layer, layer] : iter::enumerate(unitary_layers)) {
        // Generate
        if(layer.empty()) continue;
        size_t gate_size     = layer.front().pos.size();
        size_t pos_max       = std::max(layer.front().pos.back(), layer.back().pos.back());
        auto   gate_sequence = qm::get_gate_sequence(layer);
        if(deb and net.empty()) {
            auto str_rho = fmt::format("[{}]", pos_rho);
            empty_layer  = fmt::format("{0:^{1}}", " ", tw + pos_max * (uw - op) + op);
            net.emplace_back(empty_layer);
            net.back().replace(0, tw, "rho  :");
            net.back().replace(tw + pos_rho * (uw - op), str_rho.size(), str_rho);
            log.emplace_back(fmt::format("insert rho[{}] now{}", pos_rho, g.pos));
        }
        for(const auto &[idx_sublayer, seq] : iter::enumerate(gate_sequence)) {
            std::string layer_str = empty_layer;
            std::string story_str;
            for(const auto &[idx_seq, pos_gate] : iter::enumerate(seq)) {
                auto &u            = layer.at(pos_gate);
                auto  idx_sublayer = num::mod<size_t>(pos_gate, gate_size);
                if constexpr(settings::debug_circuit) layer_str.replace(0, tw, fmt::format("u[{:^2}]:", 2 * idx_layer + idx_sublayer)); // Setup layer tag

                // Going through the sequence first forward, then backward, we are handed u gates which may or may not connect to our current g gate.
                // For a successful connection, at least one pos in g should be present in u gate.
                // After connecting a full layer, trace away legs of g that are outside of the light-cone intersection between rho and sigma.

                // Check if g.pos and u.pos have sites in common
                std::vector<size_t> pos_isect;
                std::set_intersection(g.pos.begin(), g.pos.end(), u.pos.begin(), u.pos.end(), back_inserter(pos_isect));

                if(not pos_isect.empty()) {
                    // Found a matching u. Connect it
                    auto pos_old = g.pos;
                    g            = g.insert(u);
                    if constexpr(settings::debug_circuit)
                        layer_str.replace(tw + u.pos.front() * (uw - op), uw,
                                          fmt::format("[{1:^{0}}]", uw - 2, fmt::format("{:<2},{:>2}", u.pos.front(), u.pos.back())));
                    if constexpr(settings::debug_circuit) story_str.append(fmt::format("insert u{} ", u.pos));
                }
            }
            // Determine the positions that are allowed
            const std::vector<size_t> &pos_needed = lc2[2 * idx_layer + idx_sublayer + 1]; // This specifies sites that are needed to connect the coming gate

            // Check if g.pos has non-needed sites
            std::vector<size_t> pos_outside;
            std::set_difference(g.pos.begin(), g.pos.end(), pos_needed.begin(), pos_needed.end(), back_inserter(pos_outside));
            if(not pos_outside.empty()) {
                // Found positions outside of the light cone. Trace them
                auto pos_old = g.pos;
                g            = g.trace_pos(pos_outside);
                // Normalize by divinding the trace of each 2x2 identity.
                // When out of the light cone, gates do not connect to the final operator, and on those sites we
                // the unitary circuit cancels, becoming equivalent to an identity. (as in U^dagger U = I).
                for([[maybe_unused]] const auto &p : pos_outside) g.op = g.op * g.op.constant(0.5);
                if constexpr(settings::debug_circuit) story_str.append(fmt::format("trace{} ", pos_outside));
            }
            if constexpr(settings::debug_circuit) story_str.append(fmt::format("now{} ", g.pos));
            if constexpr(settings::debug_circuit) net.emplace_back(layer_str);
            if constexpr(settings::debug_circuit) log.emplace_back(story_str);
        }
    }
    cplx result;
    if(g.pos.empty()) {
        if(g.op.dimension(0) * g.op.dimension(1) != 1) // g.op should be a rank-2 tensor of dimensions 1x1
            throw except::runtime_error("Expected empty gate to have cplx op: Got dims {}", g.op.dimensions());
        // This happens when the light cone from site i can't reach site j within the current depth of the unitary circuit.
        // Essentially, ρ'_i has no support where σ^z_j connects. Therefore we effectively get
        //      Tr(ρ'_i σ^z_j) =  Tr(ρ_i ⊗ σ^z_j) =  Tr(ρ_i) Tr(σ^z_j) = 0

        Eigen::Tensor<cplx, 0> g_sig_trace = g.op.trace() * sig_gate.op.trace();
        result                             = g_sig_trace.coeff(0);
        if constexpr(settings::debug_circuit) net.emplace_back(empty_layer);
        if constexpr(settings::debug_circuit) net.back().replace(0, tw, "sig  :");
        if constexpr(settings::debug_circuit) log.emplace_back(fmt::format("sigma not connected -> result = {:.1f}{:+.1f}i", result.real(), result.imag()));
    } else {
        // In the last step we connect the sigma operator and trace everything down to a cplx
        auto g_trace     = g.trace();
        auto norm        = (rho * sig).trace();               // Will be 1 or 2 depending on what rho is (i.e. either sz or 0.5*(1+sz))
        auto g_sig_trace = g.connect_under(sig_gate).trace(); //
        if(std::abs(g_trace) > 1e-8) norm *= g_trace;
        result = g_sig_trace / norm; // What is the correct liom normalization here?
        if constexpr(settings::debug_circuit) {
            auto str_sig = fmt::format("[{}]", pos_sig);
            net.emplace_back(empty_layer);
            net.back().replace(0, tw, "sig  :");
            net.back().replace(tw + pos_sig * (uw - op), str_sig.size(), str_sig);
            log.emplace_back(fmt::format("insert sigma[{0}] now{1} trace{1} result = {2:.1f}{3:+.1f}i", pos_sig, g.pos, result.real(), result.imag()));
        }
    }

    tools::log->debug("Computed Tr(ρ_{} σ_{}) / Tr(ρ_{}) = {:.6f}{:+.6f}i", pos_rho, pos_sig, pos_rho, result.real(), result.imag());
    if constexpr(settings::debug_circuit)
        for(const auto &[idx, layer] : iter::enumerate_reverse(net)) tools::log->debug("{} | log: {}", layer, log[idx]);

    return result;
}

Eigen::Tensor<qm::cplx, 2> qm::lbit::get_lbit_real_overlap(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites) {
    auto                       ssites = static_cast<long>(sites);
    Eigen::Tensor<qm::cplx, 2> lbit_overlap(ssites, ssites);

    // We calculate the operator overlap = Tr(ρ' σ^z_j) / Tr(ρ') = 1/2^L Tr(σ' σ^z_j)
    // Define
    //      ρ = (1/2^L)  (σ^z_i + 1 ) ⊗ I_i',
    // where I_i' is the identity matrix on sites j != i
    // Then ρ represents a state |psi><psi| where site i is at level 0 with probability 1, and
    // the others sites are at level 0 or 1 with probability 1/2.
    // Eg. ρ could be a state with magnetization 1 at site i, and 0 elsewhere,
    // or ρ could be a state with 1 particle at site i, and cat state of 0 and 1 particles elsewhere.
    Eigen::MatrixXcd rho = 0.5 * (qm::spin::half::sz + qm::spin::half::id);
#pragma omp parallel for collapse(2) schedule(dynamic)
    for(long j = 0; j < ssites; j++) {
        for(long i = 0; i < ssites; i++) {
            lbit_overlap(i, j) = qm::lbit::get_lbit_exp_value(unitary_layers, rho, static_cast<size_t>(i), qm::spin::half::sz, static_cast<size_t>(j));
        }
    }
    // We require that lbit_overlap(i,j) has rows that sum up to 1
    auto sums = tenx::MatrixMap(lbit_overlap).rowwise().sum();
    if(not sums.cwiseAbs().isOnes(1e-4)) {
        tools::log->info("lbit_overlap: \n{}\nsums\n{}\n", linalg::tensor::to_string(lbit_overlap, 6), linalg::matrix::to_string(sums, 6));
        throw except::logic_error("lbit overlap rows do not sum to one. Perhaps normalization is wrong");
    }

    return lbit_overlap;
}

Eigen::Tensor<cplx, 2> qm::lbit::get_lbit_overlap_permuted(const Eigen::Tensor<cplx, 2> &lbit_overlap) {
    // First, subtract the center position of each lbit, so we get L lbits centered around zero.
    // In practice, we make a cyclic permutation of the rows of lbit_overlap
    // In addition, we mirror the lbit along its vertical, so that we can average its left and right half together
    long                   rows = lbit_overlap.dimension(0);
    long                   cols = lbit_overlap.dimension(1);
    Eigen::Tensor<cplx, 2> lbit_overlap_permuted(rows, cols);
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

std::tuple<double, double, std::vector<double>, size_t> qm::lbit::get_characteristic_length_scale(const Eigen::Tensor<cplx, 2> &lbit_overlap_permuted) {
    // Average along each column to get an estimate of the lbit
    Eigen::Tensor<double, 1> lbit_overlap_avg = lbit_overlap_permuted.real().mean(std::array<long, 1>{0});
    Eigen::Tensor<double, 1> lbit_overlap_log = lbit_overlap_avg.log();
    // Data becomes noisy if the exponential has decayed, so find a cutoff to get the slope using only the first part of the curve
    auto yavg         = std::vector<double>(lbit_overlap_avg.data(), lbit_overlap_avg.data() + lbit_overlap_avg.size());
    auto v            = stat::find_last_valid_point(yavg);
    auto c            = std::count_if(yavg.begin(), yavg.begin() + static_cast<long>(v), [](auto &val) { return val > 1e-16; });
    auto x            = num::range<double>(0, c);
    auto ylog         = std::vector<double>(lbit_overlap_log.data(), lbit_overlap_log.data() + c);
    auto [slope, res] = stat::slope(x, ylog, 0, c);
    double cls        = 1.0 / std::abs(slope);
    tools::log->debug("Computed lbit cls {:>8.6f} | sse {:>8.6f} | using {} points", cls, res, c);
    return {cls, res, yavg, c};
}

Eigen::Tensor<cplx, 2> qm::lbit::get_lbit_overlap_averaged(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_overlap_vec) {
    Eigen::Tensor<cplx, 2> avg;
    //    Eigen::Tensor<cplx,2> err;
    if(not lbit_overlap_vec.empty()) {
        long              rows = lbit_overlap_vec.front().dimension(0);
        long              cols = lbit_overlap_vec.front().dimension(1);
        size_t            reps = lbit_overlap_vec.size();
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

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::Tensor<double, 3>, Eigen::Tensor<double, 4>>
    qm::lbit::get_lbit_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec, size_t sites, size_t reps) {
    auto t_lbit_analysis = tid::tic_scope("lbit_analysis");

    long                     rows = static_cast<long>(fmix_vec.size());
    long                     cols = static_cast<long>(udepth_vec.size());
    Eigen::MatrixXd          cls_avg(rows, cols);
    Eigen::MatrixXd          cls_err(rows, cols);
    Eigen::MatrixXd          sse_avg(rows, cols);
    Eigen::MatrixXd          sse_err(rows, cols);
    Eigen::Tensor<double, 3> lbit_decay(rows, cols, static_cast<long>(sites));
    Eigen::Tensor<double, 4> lbit_lioms(rows, cols, static_cast<long>(sites), static_cast<long>(sites));
    lbit_decay.setZero();
    lbit_lioms.setZero();
    std::array<long, 3> offset3{}, extent3{};
    std::array<long, 4> offset4{}, extent4{};

    for(size_t uidx = 0; uidx < udepth_vec.size(); uidx++) {
        for(size_t fidx = 0; fidx < fmix_vec.size(); fidx++) {
            //    for(const auto & [uidx,udep] : iter::enumerate(udepth_vec)){

            //        for (const auto & [fidx,fmix] : iter::enumerate(fmix_vec) ){
            auto ur_iter = tid::ur("lbit-cls");
            ur_iter.tic();
            auto                                fmix = fmix_vec[fidx];
            auto                                udep = udepth_vec[uidx];
            std::vector<double>                 cls_vec(reps);
            std::vector<double>                 sse_vec(reps);
            std::vector<Eigen::Tensor<cplx, 2>> lbit_overlap_vec(reps);
            for(size_t i = 0; i < reps; i++) {
                std::vector<std::vector<qm::Gate>> layers;
                for(size_t l = 0; l < udep; l++) layers.emplace_back(qm::lbit::get_unitary_2gate_layer(sites, fmix));
                lbit_overlap_vec[i] = qm::lbit::get_lbit_real_overlap(layers, sites);
            }

            auto lbit_overlap_avg = qm::lbit::get_lbit_overlap_averaged(lbit_overlap_vec);
            auto lbit_overlap_per = qm::lbit::get_lbit_overlap_permuted(lbit_overlap_avg);

            auto [cls, sse, y, c] = qm::lbit::get_characteristic_length_scale(lbit_overlap_per);
            ur_iter.toc();

#if defined(_OPENMP)
            tools::log->info("Computed u {} | f {:.4f} | lbit cls {:>8.6f} | sse {:>8.6f} | threads {} | time {:8.3f} s | decay {:2} sites: {:8.2e}", udep,
                             fmix, cls, sse, omp_get_num_threads(), ur_iter.get_last_interval(), c, fmt::join(y, ", "));
#else
            tools::log->info("Computed u {} | f {:.4f} | lbit cls {:>8.6f} | sse {:>8.6f} | time {:8.3f} s | decay {:2} sites: {:8.2e}", udep, fmix, cls, sse,
                             ur_iter.get_last_interval(), c, fmt::join(y, ", "));
#endif
            cls_avg(static_cast<long>(fidx), static_cast<long>(uidx)) = cls;
            sse_avg(static_cast<long>(fidx), static_cast<long>(uidx)) = sse;

            offset3                            = {static_cast<long>(fidx), static_cast<long>(uidx), 0};
            extent3                            = {1, 1, static_cast<long>(y.size())};
            offset4                            = {static_cast<long>(fidx), static_cast<long>(uidx), 0, 0};
            extent4                            = {1, 1, lbit_overlap_avg.dimension(0), lbit_overlap_avg.dimension(1)};
            lbit_decay.slice(offset3, extent3) = Eigen::TensorMap<Eigen::Tensor<double, 3>>(y.data(), extent3);
            lbit_lioms.slice(offset4, extent4) = lbit_overlap_avg.real().reshape(extent4);
        }
    }
    return {cls_avg, sse_avg, lbit_decay, lbit_lioms};

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
