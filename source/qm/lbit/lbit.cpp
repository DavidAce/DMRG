#include "../lbit.h"
#include "../spin.h"
#include "math/tenx.h"
#include <algorithm>
#include <config/debug.h>
#include <general/iter.h>
#include <io/fmt.h>
#include <io/spdlog.h>
#include <math/linalg.h>
#include <math/num.h>
#include <math/rnd.h>
#include <math/stat.h>
#include <tools/common/log.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

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
        double               th0 = rnd::uniform_double_box(1, -1);
        double               th1 = rnd::uniform_double_box(1, -1);
        double               th2 = rnd::uniform_double_box(1, -1);
        double               th3 = rnd::uniform_double_box(1, -1);
        std::complex<double> t(rnd::uniform_double_box(1, -1), rnd::uniform_double_box(1, -1));
        auto                 indices = std::vector<size_t>{idx, idx + 1};
        Eigen::Matrix4cd     H       = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
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
            //         0                   1      0              0      1               0
            //         |                   |      |              |      |               |
            //   [ exp(-ifH) ]  --->    [ exp(-ifH) ]   --->  [ exp(-ifH) ]  --->  [ exp(-ifH) ]
            //        |                   |      |              |      |                |
            //        1                   3      2              2      3                1
            Eigen::Tensor<cplx, 2> H_shuffled = tenx::TensorMap(H, 2, 2, 2, 2).shuffle(tenx::array4{1, 0, 3, 2}).reshape(tenx::array2{4, 4});
            Eigen::MatrixXcd       expifH     = (imn * fmix * tenx::MatrixMap(H_shuffled)).exp();
            unitaries.emplace_back(tenx::TensorMap(expifH), indices, spin_dims);
        }
    }
    if constexpr(settings::debug){
        // Sanity check
        for(const auto &u : unitaries)
            if(not tenx::MatrixMap(u.op).isUnitary()) throw std::logic_error("u is not unitary!");
    }

    return unitaries;
}

Eigen::Tensor<cplx, 2> qm::lbit::get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<cplx, 2> &H) {
    // Given a matrix H, this returns exp(delta_t * H)
    // For time evolution, just make sure delta_t = -i*d,  where d is a (small) real positive number.
    return tenx::TensorCast((delta_t * tenx::MatrixMap(H)).exp());
}

std::vector<qm::Gate> qm::lbit::get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite) {
    std::vector<Gate> time_evolution_gates;
    time_evolution_gates.reserve(hams_nsite.size());
    for(auto &h : hams_nsite) time_evolution_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
    if constexpr (settings::debug){
        for(auto &t : time_evolution_gates)
            if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                throw std::runtime_error(fmt::format("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op)));
            }
    }

    return time_evolution_gates;
}

std::vector<Eigen::Tensor<cplx, 2>> qm::lbit::get_time_evolution_operators_2site(size_t sites, cplx delta_t,
                                                                                 const std::vector<Eigen::Tensor<cplx, 2>> &hams_2site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{j,j+1} and h_{j+1,j+2} commute. Therefore we can immediately use the relation
    //      exp(-i*dt *[h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}]) =  exp(-i*dt [h_{j,j+1}]) * exp(-i*dt*[h_{j+1,j+2}]) * ... * exp(-i*dt*[h_{L-2, L-1}])
    // without passing through the Suzuki-Trotter decomposition.
    // Here we expect "hams_2site" to contain terms like  h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}.

    if(hams_2site.size() != sites - 1)
        throw std::logic_error(fmt::format("Wrong number of twosite hamiltonians: {}. Expected {}", hams_2site.size(), sites - 1));

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

    if(hams_3site.size() != sites - 2)
        throw std::logic_error(fmt::format("Wrong number of three-site hamiltonians: {}. Expected {}", hams_3site.size(), sites - 2));

    std::vector<Eigen::Tensor<cplx, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_3site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

qm::cplx qm::lbit::get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &tau, size_t pos_tau,
                                      const Eigen::Matrix2cd &sig, size_t pos_sig) {
    // Generate gates for the operators
    tools::log->trace("Computing Trace (tau_{} sig_{})", pos_tau, pos_sig);
    //    auto tau_gate = qm::Gate{tau, {pos_tau}, {2l}};
    //    auto sig_gate = qm::Gate{sig, {pos_sig}, {2l}};
    auto                     tau_gate = qm::Gate(tau, {pos_tau}, {2l});
    auto                     sig_gate = qm::Gate(sig, {pos_sig}, {2l});
    auto                     g        = tau_gate; // Start with the bottom tau gate
    auto                     lc2      = qm::get_lightcone_intersection(unitary_layers, pos_tau, pos_sig);
    bool                     deb      = tools::log->level() <= spdlog::level::debug;
    std::vector<std::string> net; // Great for debugging
    std::vector<std::string> log; // Great for debugging
    std::string              empty_layer;
    size_t                   uw = 7; // Width of a unitary 2-site gate box
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
    cplx result;
    if(g.pos.empty()) {
        if(g.op.dimension(0) * g.op.dimension(1) != 1)
            throw std::runtime_error(fmt::format("Expected empty gate to have cplx op: Got dims {}", g.op.dimensions()));
        result = g.op.coeff(0);
        if(deb) net.emplace_back(empty_layer);
        if(deb) net.back().replace(0, tw, "sig  :");
        if(deb) log.emplace_back(fmt::format("sigma not connected -> result = {:.1f}{:+.1f}i", result.real(), result.imag()));
    } else {
        // In the last step we connect the sigma operator and trace everything down to a cplx
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

Eigen::Tensor<qm::cplx, 2> qm::lbit::get_lbit_real_overlap(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites) {
    Eigen::Tensor<qm::cplx, 2> lbit_overlap;
    lbit_overlap.resize(static_cast<long>(sites), static_cast<long>(sites));
    for(long j = 0; j < lbit_overlap.dimension(1); j++)
        for(long i = 0; i < lbit_overlap.dimension(0); i++)
            lbit_overlap(i, j) =
                qm::lbit::get_lbit_exp_value(unitary_layers, qm::spin::half::sz, static_cast<size_t>(i), qm::spin::half::sz, static_cast<size_t>(j));
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

#pragma omp parallel for collapse(2) schedule(dynamic)
    for(size_t uidx = 0; uidx < udepth_vec.size(); uidx++) {
        for(size_t fidx = 0; fidx < fmix_vec.size(); fidx++) {
            //    for(const auto & [uidx,udep] : iter::enumerate(udepth_vec)){

            //        for (const auto & [fidx,fmix] : iter::enumerate(fmix_vec) ){
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
#if defined(_OPENMP)
            tools::log->info("Computed u {} | f {:.4f} | lbit width {:.6f} | sse {:.6f} | threads {} | points {}: {:.8f}", udep, fmix, cls, sse,
                             omp_get_num_threads(), c, fmt::join(y, ", "));
#else
            tools::log->info("Computed u {} | f {:.4f} | lbit width {:.6f} | sse {:.6f} | points {}: {:.8f}", udep, fmix, cls, sse, c, fmt::join(y, ", "));
#endif
            cls_avg(static_cast<long>(fidx), static_cast<long>(uidx)) = cls;
            sse_avg(static_cast<long>(fidx), static_cast<long>(uidx)) = sse;

            offset3                            = {static_cast<long>(fidx), static_cast<long>(uidx), 0};
            extent3                            = {1, 1, static_cast<long>(y.size())};
            offset4                            = {static_cast<long>(fidx), static_cast<long>(uidx), 0, 0};
            extent4                            = {1, 1, lbit_overlap_avg.dimension(0), lbit_overlap_avg.dimension(1)};
            lbit_decay.slice(offset3, extent3) = Eigen::TensorMap<Eigen::Tensor<double, 1>>(y.data(), y.size());
            lbit_lioms.slice(offset4, extent4) = lbit_overlap_avg.real().reshape(extent4);
        }
    }
    return {cls_avg, sse_avg, lbit_decay, lbit_lioms};
}
