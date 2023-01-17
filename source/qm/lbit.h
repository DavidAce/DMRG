#pragma once
#include "gate.h"
#include "math/tenx/fwd_decl.h"
#include "qm.h"
#include <complex>
#include <vector>

enum class UnitaryGateType;

namespace qm::lbit {
    struct UnitaryGateProperties {
        UnitaryGateType     ugate_type;
        std::vector<double> fields;
        double              fieldvar;
    };
    struct lbitSupportAnalysis {
        Eigen::MatrixXd          cls_avg; // Characteristic length-scale of lbits
        Eigen::MatrixXd          cls_err; // Standard error of cls
        Eigen::MatrixXd          sse_avg; // Squared sum error of fits
        Eigen::MatrixXd          sse_err;
        Eigen::Tensor<double, 3> decay; // The disorder-averaged and permuted decay of lbits
        Eigen::Tensor<double, 5> supps; // The raw data from l-bit support matrices O(i,j) for each realization
        Eigen::Tensor<double, 5> pupps; // The permuted lbit support matrices O(i, |i-j|) for each realization
        lbitSupportAnalysis() {
            cls_avg.setZero();
            cls_err.setZero();
            sse_avg.setZero();
            sse_err.setZero();
            decay.setZero();
            supps.setZero();
            pupps.setZero();
        }
        template<typename T1, typename T2>
        lbitSupportAnalysis(T1 rows, T1 cols, T2 reps, T2 size) : lbitSupportAnalysis() {
            static_assert(std::is_integral_v<T1>);
            static_assert(std::is_integral_v<T2>);
            auto irows = static_cast<Eigen::Index>(rows);
            auto icols = static_cast<Eigen::Index>(cols);
            auto ireps = static_cast<Eigen::Index>(reps);
            auto isize = static_cast<Eigen::Index>(size);
            cls_avg.resize(irows, icols);
            cls_err.resize(irows, icols);
            sse_avg.resize(irows, icols);
            sse_err.resize(irows, icols);
            decay.resize(irows, icols, isize);
            supps.resize(irows, icols, ireps, isize, isize);
            pupps.resize(irows, icols, ireps, isize, isize);
        }
    };

    /* clang-format off */
    extern Eigen::Tensor<cplx, 2>               get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<cplx, 2> &hamiltonian);
    extern std::vector<qm::Gate>                get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite, double id_threshold = std::numeric_limits<double>::epsilon());
    extern std::vector<qm::SwapGate>            get_time_evolution_swap_gates(cplx delta_t, const std::vector<qm::SwapGate> &hams_nsite, double id_threshold = std::numeric_limits<double>::epsilon());
    extern std::vector<qm::Gate>                get_unitary_2gate_layer(size_t sites, double fmix);
    extern std::vector<qm::Gate>                get_unitary_2gate_layer_blocked(size_t sites, double fmix, const std::vector<double> & fields, double fieldvar);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_time_evolution_operators_2site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &twosite_hams);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_time_evolution_operators_3site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &hams_3site);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_time_evolution_mpos(cplx delta_t, const std::vector<Eigen::Tensor<cplx, 4>> &mpos);
    extern cplx                                 get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &rho, size_t pos_rho, const Eigen::Matrix2cd &sig, size_t pos_sig);
    extern cplx                                 get_lbit_exp_value2(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern cplx                                 get_lbit_exp_value3(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern Eigen::Tensor<cplx, 2>               get_lbit_support(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites);
    extern Eigen::Tensor<cplx, 2>               get_lbit_support_averaged(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_overlap_vec);
    extern Eigen::Tensor<cplx, 2>               get_lbit_permuted_support(const Eigen::Tensor<cplx, 2> &lbit_support);
    extern std::tuple<double,
                      double,
                      std::vector<double>,
                      size_t>                   get_characteristic_length_scale(const Eigen::Tensor<cplx, 2> &lbit_overlap_permuted);
    extern std::tuple<double,
                      double,
                      std::vector<double>,
                      size_t>                   get_characteristic_length_scale2(const Eigen::Tensor<cplx, 2> &lbit_overlap_permuted);
    extern lbitSupportAnalysis                  get_lbit_support_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec,
                                                                          size_t reps, size_t sites, const UnitaryGateProperties & ugate_props,
                                                                          bool randomize_fields = false);
    /* clang-format on */
}