#define FMT_USE_FLOAT128 1

#include "algorithms/flbit.h"
#include "config/parse.h"
#include "config/settings.h"
#include "config/threading.h"
#include "debug/exceptions.h"
#include "debug/info.h"
#include "debug/stacktrace.h"
#include "env/environment.h"
#include "io/filesystem.h"
#include "math/float.h"
#include "math/linalg.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/lbit.h"
#include "qm/spin.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/mps.h"
#include <fmt/format.h>
#include <h5pp/h5pp.h>
#include <unsupported/Eigen/MPRealSupport>

auto get_hamiltonian(const flbit &f) {
    //    for(const auto &field : f.tensors.model->get_parameter("J1_rand"))
    using matrix_t = Eigen::Matrix<cplx_t, Eigen::Dynamic, Eigen::Dynamic>;
    auto L         = f.tensors.get_length<size_t>();
    auto N         = static_cast<long>(std::pow(2.0, L));
    auto H         = matrix_t(N, N);
    auto SZ        = qm::spin::gen_manybody_spins(qm::spin::half::sz, L, true);

    auto J1 = f.tensors.model->get_parameter("J1_rand");
    auto J2 = f.tensors.model->get_parameter("J2_rand");
    auto J3 = f.tensors.model->get_parameter("J3_rand");

    H.setZero();
    for(size_t i = 0; i < L - 0; ++i) H += std::any_cast<real_t>(J1[i]) * SZ[i].cast<cplx_t>();
    for(size_t i = 0; i < L - 1; ++i) {
        auto J2i = std::any_cast<h5pp::varr_t<real_t>>(J2[i]);
        for(size_t j = i + 1; j < L; ++j) {
            auto r = j - i;
            H += J2i[r] * SZ[i].cast<cplx_t>() * SZ[j].cast<cplx_t>();
        }
    }
    for(size_t i = 0; i < L - 2; ++i) H += std::any_cast<real_t>(J3[i]) * SZ[i + 0].cast<cplx_t>() * SZ[i + 1].cast<cplx_t>() * SZ[i + 2].cast<cplx_t>();

    return H;
}

Eigen::MatrixXcd get_time_evolution_from_nbody_ham_gates(const flbit &f) {
    auto L = f.tensors.get_length<size_t>();
    auto N = static_cast<long>(std::pow(2.0, L));
    auto H = Eigen::MatrixXcd(N, N);
    auto I = Eigen::MatrixXcd(N, N);

    H.setZero();
    I.setIdentity();
    auto pos    = num::range<size_t>(0, L);
    auto dim    = std::vector<long>(L, 2);
    auto H_gate = qm::Gate(H, pos, dim);
    auto I_gate = qm::Gate(I, pos, dim);

    bool has_slow_gates = not f.ham_gates_1body.empty() or not f.ham_gates_2body.empty() or not f.ham_gates_3body.empty();
    bool has_swap_gates = not f.ham_swap_gates_1body.empty() or not f.ham_swap_gates_2body.empty() or not f.ham_swap_gates_3body.empty();
    if(has_slow_gates) {
        auto sequence1 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.ham_gates_1body, CircOp::NONE);
        auto sequence2 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.ham_gates_2body, CircOp::NONE);
        auto sequence3 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.ham_gates_3body, CircOp::NONE);
        for(const auto &idx : sequence1) {
            auto &g = f.ham_gates_1body.at(idx);
            tools::log->info("H_1body {} pos {}: {}", idx, g.pos, linalg::matrix::to_string(tenx::MatrixMap(g.op).diagonal().transpose(), 16));
            H_gate = qm::Gate(H_gate.op + g.connect_above(I_gate).op, H_gate.pos, H_gate.dim);
        }
        for(const auto &idx : sequence2) {
            auto &g = f.ham_gates_2body.at(idx);
            //            H_gate  = f.ham_gates_2body.at(idx).connect_above(H_gate);
        }
        for(const auto &idx : sequence3) {
            auto &g = f.ham_gates_3body.at(idx);
            //            H_gate  = qm::Gate(H_gate.op + g.connect_above(I_gate).op, H_gate.pos, H_gate.dim);
        }
    } else if(has_swap_gates) {
        auto sequence1 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.ham_swap_gates_1body, CircOp::NONE);
        auto sequence2 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.ham_swap_gates_2body, CircOp::NONE);
        auto sequence3 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.ham_swap_gates_3body, CircOp::NONE);
        for(const auto &idx : sequence1) {
            auto &g = f.ham_swap_gates_1body.at(idx);
            tools::log->info("H_1body {} pos {}: {}", idx, g.pos, linalg::matrix::to_string(tenx::MatrixMap(g.op).diagonal().transpose(), 16));
            H_gate = qm::Gate(H_gate.op + g.connect_above(I_gate).op, H_gate.pos, H_gate.dim);
        }
        for(const auto &idx : sequence2) {
            auto &g = f.ham_swap_gates_2body.at(idx);
            //            H_gate  = g.connect_above(H_gate);
        }
        for(const auto &idx : sequence3) {
            auto &g = f.ham_swap_gates_3body.at(idx);
            //            H_gate  = qm::Gate(g.op, g.pos, g.dim).connect_above(H_gate);
            //            H_gate = qm::Gate(H_gate.op + g.connect_above(I_gate).op, H_gate.pos, H_gate.dim);
        }
    }
    tools::log->info("H_hbody L: {}", linalg::matrix::to_string(tenx::MatrixMap(H_gate.op).diagonal().transpose(), 16));
    using namespace std::complex_literals;
    auto T_gate = H_gate.exp(cplx_t(-1.0i) * f.status.delta_t);
    return tenx::MatrixMap(T_gate.op);
}

Eigen::MatrixXcd get_time_evolution_from_nbody_time_gates(const flbit &f) {
    //    for(const auto &field : f.tensors.model->get_parameter("J1_rand"))
    auto L = f.tensors.get_length<size_t>();
    auto N = static_cast<long>(std::pow(2.0, L));
    auto T = Eigen::MatrixXcd(N, N);
    T.setIdentity();
    auto pos    = num::range<size_t>(0, L);
    auto dim    = std::vector<long>(L, 2);
    auto T_gate = qm::Gate(tenx::TensorMap(T), pos, dim);

    bool has_slow_gates = not f.time_gates_1body.empty() or not f.time_gates_2body.empty() or not f.time_gates_3body.empty();
    bool has_swap_gates = not f.time_swap_gates_1body.empty() or not f.time_swap_gates_2body.empty() or not f.time_swap_gates_3body.empty();
    if(has_slow_gates) {
        auto sequence1 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_gates_1body, CircOp::NONE);
        auto sequence2 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_gates_2body, CircOp::NONE);
        auto sequence3 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_gates_3body, CircOp::NONE);
        for(const auto &idx : sequence1) {
            auto &g = f.time_gates_1body.at(idx);
            tools::log->info("T_1body {} pos {}: {} | {}", idx, g.pos, linalg::matrix::to_string(tenx::MatrixMap(T_gate.op).diagonal().transpose(), 16),
                             linalg::matrix::to_string(tenx::MatrixMap(g.op).diagonal().transpose(), 16));

            //            T_gate    = g.connect_above(T_gate);
            if(idx == 0)
                T_gate = g;
            else {
                auto dim2 = tenx::array2{T_gate.op.dimension(0) * g.op.dimension(0), T_gate.op.dimension(1) * g.op.dimension(1)};
                auto dim  = T_gate.dim;
                auto pos  = T_gate.pos;
                for(auto p : g.pos) pos.emplace_back(p);
                for(auto d : g.dim) dim.emplace_back(d);
                T_gate = qm::Gate(T_gate.op.contract(g.op, tenx::idx()).shuffle(tenx::array4{0, 2, 1, 3}).reshape(dim2), pos, dim);
            }
        }
        //        for(const auto &idx : sequence2) T_gate = f.time_gates_2body.at(idx).connect_above(T_gate);
        //        for(const auto &idx : sequence3) T_gate = f.time_gates_3body.at(idx).connect_above(T_gate);
    } else if(has_swap_gates) {
        auto sequence1 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_swap_gates_1body, CircOp::NONE);
        auto sequence2 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_swap_gates_2body, CircOp::NONE);
        auto sequence3 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_swap_gates_3body, CircOp::NONE);
        for(const auto &idx : sequence1) {
            auto &g = f.time_swap_gates_1body.at(idx);
            T_gate  = qm::Gate(g.op, g.pos, g.dim).connect_above(T_gate);
        }
        //        for(const auto &idx : sequence2) {
        //            auto &g = f.time_swap_gates_2body.at(idx);
        //            T_gate  = g.connect_above(T_gate);
        //        }
        //        for(const auto &idx : sequence3) {
        //            auto &g = f.time_swap_gates_3body.at(idx);
        //            T_gate  = qm::Gate(g.op, g.pos, g.dim).connect_above(T_gate);
        //        }
    }
    tools::log->info("T_1body L pos L: {}", linalg::matrix::to_string(tenx::MatrixMap(T_gate.op).diagonal().transpose(), 16));
    return tenx::MatrixMap(T_gate.op);
}

Eigen::MatrixXcd get_time_evolution_from_nbody_time_matrix(const flbit &f) {
    //    for(const auto &field : f.tensors.model->get_parameter("J1_rand"))
    auto L = f.tensors.get_length<size_t>();
    auto N = static_cast<long>(std::pow(2.0, L));
    auto T = Eigen::MatrixXcd(N, N);
    auto I = Eigen::MatrixXcd(N, N);
    T.setIdentity();
    I.setIdentity();
    auto pos = num::range<size_t>(0, L);
    auto dim = std::vector<long>(L, 2);

    auto T_gate = qm::Gate(tenx::TensorMap(T), pos, dim);
    auto I_gate = qm::Gate(tenx::TensorMap(I), pos, dim);

    bool has_slow_gates = not f.time_gates_1body.empty() or not f.time_gates_2body.empty() or not f.time_gates_3body.empty();
    bool has_swap_gates = not f.time_swap_gates_1body.empty() or not f.time_swap_gates_2body.empty() or not f.time_swap_gates_3body.empty();
    if(has_slow_gates) {
        auto sequence1 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_gates_1body, CircOp::NONE);
        auto sequence2 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_gates_2body, CircOp::NONE);
        auto sequence3 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_gates_3body, CircOp::NONE);
        for(const auto &idx : sequence1) {
            auto &g = f.time_gates_1body.at(idx);
            T_gate  = qm::Gate(T_gate.op * g.connect_above(I_gate).op, T_gate.pos, T_gate.dim);
            auto m  = tenx::MatrixMap(g.op);
            tools::log->info("M_1body {} pos {}: {}", idx, g.pos, linalg::matrix::to_string(m.diagonal().transpose(), 16));
        }
        //        for(const auto &idx : sequence2) T_gate = f.time_gates_2body.at(idx).connect_above(T_gate);
        //        for(const auto &idx : sequence3) T_gate = f.time_gates_3body.at(idx).connect_above(T_gate);
    } else if(has_swap_gates) {
        auto sequence1 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_swap_gates_1body, CircOp::NONE);
        auto sequence2 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_swap_gates_2body, CircOp::NONE);
        auto sequence3 = tools::finite::mps::generate_gate_sequence(*f.state_lbit, f.time_swap_gates_3body, CircOp::NONE);
        for(const auto &idx : sequence1) {
            auto &g = f.time_swap_gates_1body.at(idx);
            T_gate  = qm::Gate(T_gate.op * g.connect_above(I_gate).op, T_gate.pos, T_gate.dim);
            auto m  = tenx::MatrixMap(g.op);
            tools::log->info("M_1body {} pos {}: {}", idx, g.pos, linalg::matrix::to_string(m.diagonal().transpose(), 16));
        }
        //        for(const auto &idx : sequence2) {
        //            auto &g = f.time_swap_gates_2body.at(idx);
        //            T_gate  = g.connect_above(T_gate);
        //        }
        //        for(const auto &idx : sequence3) {
        //            auto &g = f.time_swap_gates_3body.at(idx);
        //            T_gate  = qm::Gate(g.op, g.pos, g.dim).connect_above(T_gate);
        //        }
    }
    tools::log->info("M_1body L pos L: {}", linalg::matrix::to_string(T.diagonal().transpose(), 16));
    return tenx::MatrixMap(T_gate.op);
}

size_t assert_lbit_evolution(const flbit &f) {
    // Test time evolution as  |psi(t+dt)> = U^-1 exp(-iHdt) U |psi(t)>, starting from t = 0
    // We get a unitary tensor whose index 0 connects onto index 0 on an MPS, meaning,
    /*
     * @verbatim
     *     1---A---2
     *         |
     *         0
     *         1
     *         |
     *         U
     *         |
     *         0
     * @endverbatim
     *
     *  When we use matrix operations, the multiplication works as
     *
     *  0---U---1 0---A
     *
     *
     */
    for(const auto &mps : f.state_lbit->mps_sites) mps->assert_normalized();
    f.state_lbit->assert_validity();
    bool has_time_swap_gates = not f.time_swap_gates_Lbody.empty();
    bool has_time_slow_gates = not f.time_gates_Lbody.empty();
    bool has_unitary_gates   = not f.unitary_gates_2site_layers.empty();

    if(not has_time_swap_gates and not has_time_slow_gates) return 0;
    if(not has_unitary_gates) return 0;

    auto H_exact    = get_hamiltonian(f);
    auto H_Lbody    = not f.ham_gates_Lbody.empty() ? tenx::MatrixMap(f.ham_gates_Lbody.front().op_t) : tenx::MatrixMap(f.ham_swap_gates_Lbody.front().op_t);
    auto digits10_d = std::numeric_limits<real>::digits10;
#if __has_include(<quadmath.h>)
    auto digits10_t = FLT128_DIG;
#else
    auto digits10_t = std::numeric_limits<real_t>::digits10;
#endif
    auto epsilon_d = std::numeric_limits<real>::epsilon();
    //    auto epsilon_t  = std::numeric_limits<real_t>::epsilon();
    if(not H_exact.isApprox(H_Lbody))
        throw except::runtime_error("Hamiltonians do not match: \nH_exact: {}\nH_Lbody: {}",
                                    linalg::matrix::to_string(H_exact.diagonal().transpose().real(), digits10_t),
                                    linalg::matrix::to_string(H_Lbody.diagonal().transpose().real(), digits10_t));
    tools::log->info("H_exact: {}", linalg::matrix::to_string(H_exact.diagonal().transpose().real(), digits10_t));
    tools::log->info("H_Lbody: {}", linalg::matrix::to_string(H_Lbody.diagonal().transpose().real(), digits10_t));
    auto T_Lbody = has_time_swap_gates ? tenx::MatrixMap(f.time_swap_gates_Lbody.front().op) : tenx::MatrixMap(f.time_gates_Lbody.front().op);
    auto T_tbody = get_time_evolution_from_nbody_time_gates(f);
    auto T_mbody = get_time_evolution_from_nbody_time_matrix(f);
    auto T_hbody = get_time_evolution_from_nbody_ham_gates(f);

    tools::log->info("T_Lbody: {}", linalg::matrix::to_string(T_Lbody.diagonal().transpose(), digits10_d));
    tools::log->info("T_tbody: {}", linalg::matrix::to_string(T_tbody.diagonal().transpose(), digits10_d));
    tools::log->info("T_mbody: {}", linalg::matrix::to_string(T_mbody.diagonal().transpose(), digits10_d));
    tools::log->info("T_hbody: {}", linalg::matrix::to_string(T_hbody.diagonal().transpose(), digits10_d));
    double max_iHdt_error = 10 * epsilon_d * std::max(1.0, static_cast<double>(f.status.delta_t.real()));

    if(not T_Lbody.isApprox(T_tbody, max_iHdt_error))
        throw except::runtime_error("Time evolution operators do not match: \nT_Lbody: {}\nT_tbody: {}",
                                    linalg::matrix::to_string(T_Lbody.diagonal().transpose(), digits10_d),
                                    linalg::matrix::to_string(T_tbody.diagonal().transpose(), digits10_d));

    if(not T_Lbody.isApprox(T_hbody, max_iHdt_error))
        throw except::runtime_error("Time evolution operators do not match: \nT_Lbody: {}\nT_hbody: {}",
                                    linalg::matrix::to_string(T_Lbody.diagonal().transpose(), digits10_d),
                                    linalg::matrix::to_string(T_hbody.diagonal().transpose(), digits10_d));

    double max_time_error = 1e-2; // * f.tensors.get_length<double>() * settings::flbit::time_gate_id_threshold;
    double max_norm_error = 1e-14;
    double max_untary_err = 1e-14;

    auto unitarytensor = qm::lbit::get_unitary_circuit_as_tensor(f.unitary_gates_2site_layers);
    auto unitarymatrix = tenx::MatrixMap(unitarytensor);
    auto exp_lbit_iHdt = has_time_swap_gates ? tenx::MatrixMap(f.time_swap_gates_Lbody.front().op) : tenx::MatrixMap(f.time_gates_Lbody.front().op);

    tools::log->info("exp_lbit_iHdt: {}", linalg::matrix::to_string(exp_lbit_iHdt.diagonal().transpose(), digits10_d));

    // Check that the matrices are unitary
    if(not unitarymatrix.isUnitary()) throw except::logic_error("unitarymatrix is not unitary");
    if(not exp_lbit_iHdt.isUnitary()) throw except::logic_error("exp_lbit_iHdt is not unitary");

    auto             psi_real_init = tools::finite::measure::mps2tensor(*f.state_real_init);
    auto             psi_lbit_init = tools::finite::measure::mps2tensor(*f.state_lbit_init);
    auto             psi_real_tebd = tools::finite::measure::mps2tensor(*f.tensors.state);
    auto             psi_lbit_tebd = tools::finite::measure::mps2tensor(*f.state_lbit);
    Eigen::VectorXcd vec_real_init = tenx::VectorMap(psi_real_init);
    Eigen::VectorXcd vec_lbit_init = tenx::VectorMap(psi_lbit_init);
    Eigen::VectorXcd vec_real_tebd = tenx::VectorMap(psi_real_tebd);
    Eigen::VectorXcd vec_lbit_tebd = tenx::VectorMap(psi_lbit_tebd);
    Eigen::VectorXcd vec_lbit_levo = exp_lbit_iHdt * vec_lbit_init;                 // From lbit
    Eigen::VectorXcd vec_lbit_revo = exp_lbit_iHdt * unitarymatrix * vec_real_init; // From real
    Eigen::VectorXcd vec_lbit_devo = exp_lbit_iHdt.conjugate() * vec_lbit_tebd;     // From real
    Eigen::VectorXcd vec_real_levo = unitarymatrix.adjoint() * exp_lbit_iHdt * vec_lbit_init;
    Eigen::VectorXcd vec_real_revo = unitarymatrix.adjoint() * exp_lbit_iHdt * unitarymatrix * vec_real_init;
    Eigen::VectorXcd vec_real_back = unitarymatrix.adjoint() * vec_lbit_init; // Returns the initial lbit state to real basis

    auto overlap_lbit_levo_tebd = vec_lbit_levo.dot(vec_lbit_tebd);           // One if time evo works independently of U
    auto overlap_lbit_revo_tebd = vec_lbit_revo.dot(vec_lbit_tebd);           // One if time evo and U both work
    auto overlap_lbit_devo_init = vec_lbit_devo.dot(vec_lbit_init);           // One if time evo works independently of U
    auto overlap_real_levo_tebd = vec_real_levo.dot(vec_real_tebd);           // One if time evo and U both work
    auto overlap_real_revo_tebd = vec_real_revo.dot(vec_real_tebd);           // One if time evo and U both work
    auto overlap_real_back_init = vec_real_back.dot(vec_real_init);           // One if U works independently of time evo

    auto real_error_lbit_levo_tebd = std::abs(std::real(overlap_lbit_levo_tebd) - 1.0);
    auto real_error_lbit_revo_tebd = std::abs(std::real(overlap_lbit_revo_tebd) - 1.0);
    auto real_error_lbit_devo_init = std::abs(std::real(overlap_lbit_devo_init) - 1.0);
    auto real_error_real_levo_tebd = std::abs(std::real(overlap_real_levo_tebd) - 1.0);
    auto real_error_real_revo_tebd = std::abs(std::real(overlap_real_revo_tebd) - 1.0);
    auto real_error_real_back_init = std::abs(std::real(overlap_real_back_init) - 1.0);
    auto imag_error_lbit_levo_tebd = std::abs(std::imag(overlap_lbit_levo_tebd));
    auto imag_error_lbit_revo_tebd = std::abs(std::imag(overlap_lbit_revo_tebd));
    auto imag_error_lbit_devo_init = std::abs(std::imag(overlap_lbit_devo_init));
    auto imag_error_real_levo_tebd = std::abs(std::imag(overlap_real_levo_tebd));
    auto imag_error_real_revo_tebd = std::abs(std::imag(overlap_real_revo_tebd));
    auto imag_error_real_back_init = std::abs(std::imag(overlap_real_back_init));
    auto norm_error_lbit_levo_tebd = std::abs(std::abs(overlap_lbit_levo_tebd) - 1.0);
    auto norm_error_lbit_revo_tebd = std::abs(std::abs(overlap_lbit_revo_tebd) - 1.0);
    auto norm_error_lbit_devo_init = std::abs(std::abs(overlap_lbit_devo_init) - 1.0);
    auto norm_error_real_levo_tebd = std::abs(std::abs(overlap_real_levo_tebd) - 1.0);
    auto norm_error_real_revo_tebd = std::abs(std::abs(overlap_real_revo_tebd) - 1.0);
    auto norm_error_real_back_init = std::abs(std::abs(overlap_real_back_init) - 1.0);

    tools::log->info("vec_real_init          : {:.16f}", vec_real_init.norm());
    tools::log->info("vec_lbit_init          : {:.16f}", vec_lbit_init.norm());
    tools::log->info("vec_real_tebd          : {:.16f}", vec_real_tebd.norm());
    tools::log->info("vec_lbit_tebd          : {:.16f}", vec_lbit_tebd.norm());
    tools::log->info("vec_lbit_levo          : {:.16f}", vec_lbit_levo.norm());
    tools::log->info("vec_lbit_revo          : {:.16f}", vec_lbit_revo.norm());
    tools::log->info("vec_lbit_devo          : {:.16f}", vec_lbit_devo.norm());
    tools::log->info("vec_real_levo          : {:.16f}", vec_real_levo.norm());
    tools::log->info("vec_real_revo          : {:.16f}", vec_real_revo.norm());
    tools::log->info("vec_real_back          : {:.16f}", vec_real_back.norm());

    tools::log->info("overlap_lbit_levo_tebd : {:.16f} | {:.16f} |", overlap_lbit_levo_tebd, std::abs(overlap_lbit_levo_tebd));
    tools::log->info("overlap_lbit_revo_tebd : {:.16f} | {:.16f} |", overlap_lbit_revo_tebd, std::abs(overlap_lbit_revo_tebd));
    tools::log->info("overlap_lbit_devo_init : {:.16f} | {:.16f} |", overlap_lbit_devo_init, std::abs(overlap_lbit_devo_init));
    tools::log->info("overlap_real_levo_tebd : {:.16f} | {:.16f} |", overlap_real_levo_tebd, std::abs(overlap_real_levo_tebd));
    tools::log->info("overlap_real_revo_tebd : {:.16f} | {:.16f} |", overlap_real_revo_tebd, std::abs(overlap_real_revo_tebd));
    tools::log->info("overlap_real_lbit_init : {:.16f} | {:.16f} |", overlap_real_back_init, std::abs(overlap_real_back_init));

    /* clang-format off */
    if(std::abs(vec_real_init.norm()-1.0) > max_norm_error) throw except::logic_error("vec_real_init.norm() > {:.3e}: {:.16f}", vec_real_init.norm(), max_norm_error);
    if(std::abs(vec_lbit_init.norm()-1.0) > max_norm_error) throw except::logic_error("vec_lbit_init.norm() > {:.3e}: {:.16f}", vec_lbit_init.norm(), max_norm_error);
    if(std::abs(vec_real_tebd.norm()-1.0) > max_norm_error) throw except::logic_error("vec_real_tebd.norm() > {:.3e}: {:.16f}", vec_real_tebd.norm(), max_norm_error);
    if(std::abs(vec_lbit_tebd.norm()-1.0) > max_norm_error) throw except::logic_error("vec_lbit_tebd.norm() > {:.3e}: {:.16f}", vec_lbit_tebd.norm(), max_norm_error);
    if(std::abs(vec_lbit_levo.norm()-1.0) > max_norm_error) throw except::logic_error("vec_lbit_levo.norm() > {:.3e}: {:.16f}", vec_lbit_levo.norm(), max_norm_error);
    if(std::abs(vec_lbit_revo.norm()-1.0) > max_norm_error) throw except::logic_error("vec_lbit_revo.norm() > {:.3e}: {:.16f}", vec_lbit_revo.norm(), max_norm_error);
    if(std::abs(vec_lbit_devo.norm()-1.0) > max_norm_error) throw except::logic_error("vec_lbit_devo.norm() > {:.3e}: {:.16f}", vec_lbit_devo.norm(), max_norm_error);
    if(std::abs(vec_real_levo.norm()-1.0) > max_norm_error) throw except::logic_error("vec_real_levo.norm() > {:.3e}: {:.16f}", vec_real_levo.norm(), max_norm_error);
    if(std::abs(vec_real_revo.norm()-1.0) > max_norm_error) throw except::logic_error("vec_real_revo.norm() > {:.3e}: {:.16f}", vec_real_revo.norm(), max_norm_error);
    if(std::abs(vec_real_back.norm()-1.0) > max_norm_error) throw except::logic_error("vec_real_back.norm() > {:.3e}: {:.16f}", vec_real_back.norm(), max_norm_error);

    if(real_error_lbit_levo_tebd > max_time_error) throw except::logic_error("real_error_lbit_levo_tebd > {:.3e}: {:.16f}", real_error_lbit_levo_tebd, max_time_error);
    if(real_error_lbit_revo_tebd > max_time_error) throw except::logic_error("real_error_lbit_revo_tebd > {:.3e}: {:.16f}", real_error_lbit_revo_tebd, max_time_error);
    if(real_error_lbit_devo_init > max_time_error) throw except::logic_error("real_error_lbit_devo_init > {:.3e}: {:.16f}", real_error_lbit_devo_init, max_time_error);
    if(real_error_real_levo_tebd > max_time_error) throw except::logic_error("real_error_real_levo_tebd > {:.3e}: {:.16f}", real_error_real_levo_tebd, max_time_error);
    if(real_error_real_revo_tebd > max_time_error) throw except::logic_error("real_error_real_revo_tebd > {:.3e}: {:.16f}", real_error_real_revo_tebd, max_time_error);
    if(real_error_real_back_init > max_time_error) throw except::logic_error("real_error_real_back_init > {:.3e}: {:.16f}", real_error_real_back_init, max_time_error);
    if(imag_error_lbit_levo_tebd > max_time_error) throw except::logic_error("imag_error_lbit_levo_tebd > {:.3e}: {:.16f}", imag_error_lbit_levo_tebd, max_time_error);
    if(imag_error_lbit_revo_tebd > max_time_error) throw except::logic_error("imag_error_lbit_revo_tebd > {:.3e}: {:.16f}", imag_error_lbit_revo_tebd, max_time_error);
    if(imag_error_lbit_devo_init > max_time_error) throw except::logic_error("imag_error_lbit_devo_init > {:.3e}: {:.16f}", imag_error_lbit_devo_init, max_time_error);
    if(imag_error_real_levo_tebd > max_time_error) throw except::logic_error("imag_error_real_levo_tebd > {:.3e}: {:.16f}", imag_error_real_levo_tebd, max_time_error);
    if(imag_error_real_revo_tebd > max_time_error) throw except::logic_error("imag_error_real_revo_tebd > {:.3e}: {:.16f}", imag_error_real_revo_tebd, max_time_error);
    if(imag_error_real_back_init > max_time_error) throw except::logic_error("imag_error_real_back_init > {:.3e}: {:.16f}", imag_error_real_back_init, max_time_error);
    if(norm_error_lbit_levo_tebd > max_time_error) throw except::logic_error("norm_error_lbit_levo_tebd > {:.3e}: {:.16f}", norm_error_lbit_levo_tebd, max_time_error);
    if(norm_error_lbit_revo_tebd > max_time_error) throw except::logic_error("norm_error_lbit_revo_tebd > {:.3e}: {:.16f}", norm_error_lbit_revo_tebd, max_time_error);
    if(norm_error_lbit_devo_init > max_time_error) throw except::logic_error("norm_error_lbit_devo_init > {:.3e}: {:.16f}", norm_error_lbit_devo_init, max_time_error);
    if(norm_error_real_levo_tebd > max_time_error) throw except::logic_error("norm_error_real_levo_tebd > {:.3e}: {:.16f}", norm_error_real_levo_tebd, max_time_error);
    if(norm_error_real_revo_tebd > max_time_error) throw except::logic_error("norm_error_real_revo_tebd > {:.3e}: {:.16f}", norm_error_real_revo_tebd, max_time_error);
    if(norm_error_real_back_init > max_untary_err) throw except::logic_error("norm_error_real_back_init > {:.3e}: {:.16f}", norm_error_real_back_init, max_untary_err);
    /* clang-format on */

    return 1;
}

void set_log(const std::string &name) { tools::log = tools::Logger::setLogger(name, settings::console::loglevel, settings::console::timestamp); }

int main(int argc, char *argv[]) {
    settings::parse(argc, argv);

    tools::log = tools::Logger::setLogger("lbit-evoution-test", settings::console::loglevel, settings::console::timestamp);

    // Set up the number of openmp and std threads for Eigen Tensor
    settings::configure_threads();

    // Seed with random::device initially (This also takes care of srand used by Eigen)
    // This is to make reproducible simulations
    rnd::seed(settings::input::seed);

    // Register termination codes and what to do in those cases
    debug::register_callbacks();

    tools::log->info("real   epsilon: {:.4e}", std::numeric_limits<double>::epsilon());
    tools::log->info("real   max_dig: {}", std::numeric_limits<double>::max_digits10);
    tools::log->info("real_t epsilon: {:.4e}", std::numeric_limits<long double>::epsilon());
    tools::log->info("real_t max_dig: {}", std::numeric_limits<long double>::max_digits10);
#if __has_include(<quadmath.h>)
    tools::log->info(" __float128 max_dig: {}", FLT128_DIG);
#endif
    // Initialize the flbit algorithm
    auto flbit_swap = flbit(nullptr);            // Will evolve with swap gates on
    auto flbit_slow = flbit(nullptr);            // Will evolve with swap gates off

    settings::flbit::time_gate_id_threshold = 0; // Turn off time gate skips
    settings::flbit::use_swap_gates         = true;
    set_log("swap");
    flbit_swap.run_preprocessing();

    settings::flbit::use_swap_gates = false;
    set_log("slow");
    flbit_slow.run_preprocessing();

    size_t test_swap_count = 0;
    size_t test_slow_count = 0;

    while(flbit_swap.status.algo_stop == AlgorithmStop::NONE and flbit_slow.status.algo_stop == AlgorithmStop::NONE) {
        set_log("swap");
        flbit_swap.update_state();
        tools::finite::measure::do_all_measurements(*flbit_swap.tensors.state);
        flbit_swap.print_status();

        set_log("slow");
        flbit_slow.update_state();
        tools::finite::measure::do_all_measurements(*flbit_slow.tensors.state);
        flbit_slow.print_status();

        set_log("slow");
        test_slow_count += assert_lbit_evolution(flbit_slow);

        set_log("swap");
        test_swap_count += assert_lbit_evolution(flbit_swap);

        set_log("swap");
        flbit_swap.update_time_evolution_gates();
        set_log("slow");
        flbit_slow.update_time_evolution_gates();
    }
    if(test_swap_count == 0) { throw except::logic_error("No swap tests ran"); }
    if(test_slow_count == 0) { throw except::logic_error("No slow tests ran"); }
    if(test_swap_count != test_slow_count) { throw except::logic_error("Unequal number of tests ran: {} != {}", test_swap_count, test_slow_count); }
    return 0;
}
