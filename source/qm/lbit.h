#pragma once
#include "gate.h"
#include "math/tenx/fwd_decl.h"
#include "qm.h"
#include <complex>
#include <vector>

enum class UnitaryGateWeight;

namespace qm::lbit {
    struct UnitaryGateProperties {
        size_t                      sites; /*!< Width of the circuit, i.e. system length/number of sites */
        size_t                      depth; /*!< Number of layers in the unitary circuit (1 layer connects all neighbors once) */
        double                      fmix;  /*!< The mixing factor f in exp(-ifM) */
        double                      tstd;  /*!< Standard deviation for the theta parameters in the unitary gate */
        double                      cstd;  /*!< Standard deviation for the c parameters in the unitary gate  */
        UnitaryGateWeight           tgw8;  /*!< Choose IDENTITY|EXPDECAY type of gate weight (hvals is required for EXPDECAY) */
        UnitaryGateWeight           cgw8;  /*!< Choose IDENTITY|EXPDECAY type of gate weight (hvals is required for EXPDECAY) */
        double                      hmean; /*!< mean of random onsite fields */
        double                      hwdth; /*!< width of random onsite fields (st.dev. if normal) */
        std::string_view            hdist; /*!< distribution of onsite fields */
        mutable std::vector<double> hvals; /*!< onsite fields of the l-bit hamiltonian, needed for type == UnitaryGateWeight::EXPDECAY */
        UnitaryGateProperties() = default;
        UnitaryGateProperties(const std::vector<double> &h = {});
        std::string string() const;
        void        randomize_hvals() const;
    };
    struct lbitSupportAnalysis {
        Eigen::Tensor<double, 6> cls_avg; // Characteristic length-scale of lbits
        Eigen::Tensor<double, 6> cls_err; // Standard error of cls
        Eigen::Tensor<double, 6> sse_avg; // Squared sum error of fits
        Eigen::Tensor<double, 6> sse_err;
        Eigen::Tensor<double, 7> decay_avg; // The decay of l-bits: permuted and averaged over site and disorder
        Eigen::Tensor<double, 7> decay_err; // The sterr of l-bits: permuted and averaged over site and disorder
        Eigen::Tensor<double, 9> support;   // The raw data from l-bit support matrices O(i,j) for each realization
        Eigen::Tensor<double, 9> permute;   // The permuted lbit support matrices O(i, |i-j|) for each realization
        lbitSupportAnalysis() {
            cls_avg.setZero();
            cls_err.setZero();
            sse_avg.setZero();
            sse_err.setZero();
            decay_avg.setZero();
            decay_err.setZero();
            support.setZero();
            permute.setZero();
        }
        lbitSupportAnalysis(size_t ndpth, size_t nfmix, size_t ntstd, size_t ncstd, size_t ntgw8, size_t ncgw8, size_t nreps, size_t nsize)
            : lbitSupportAnalysis() {
            auto idpth = static_cast<Eigen::Index>(ndpth);
            auto ifmix = static_cast<Eigen::Index>(nfmix);
            auto itstd = static_cast<Eigen::Index>(ntstd);
            auto icstd = static_cast<Eigen::Index>(ncstd);
            auto itgw8 = static_cast<Eigen::Index>(ntgw8);
            auto icgw8 = static_cast<Eigen::Index>(ncgw8);
            auto ireps = static_cast<Eigen::Index>(nreps);
            auto isize = static_cast<Eigen::Index>(nsize);
            cls_avg.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8);
            cls_err.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8);
            sse_avg.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8);
            sse_err.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8);
            decay_avg.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8, isize);
            decay_err.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8, isize);
            support.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8, ireps, isize, isize);
            permute.resize(idpth, ifmix, itstd, icstd, itgw8, icgw8, ireps, isize, isize);
        }
    };

    /* clang-format off */
    extern Eigen::Tensor<cplx, 2>               get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<cplx, 2> &hamiltonian);
    extern std::vector<qm::Gate>                get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite, double id_threshold = std::numeric_limits<double>::epsilon());
    extern std::vector<qm::SwapGate>            get_time_evolution_swap_gates(cplx delta_t, const std::vector<qm::SwapGate> &hams_nsite, double id_threshold = std::numeric_limits<double>::epsilon());
//    extern std::vector<qm::Gate>                get_unitary_2gate_layer(size_t sites, double fmix);
    extern std::vector<qm::Gate>                get_unitary_2gate_layer(const UnitaryGateProperties & u);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_unitary_mpo_layer(const std::vector<qm::Gate> & ulayer);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_unitary_mpo_layer(const UnitaryGateProperties & u);
    extern std::vector<Eigen::Tensor<cplx, 4>>  merge_unitary_mpo_layers(const std::vector<Eigen::Tensor<cplx, 4>> & mpos_dn, const std::vector<Eigen::Tensor<cplx, 4>> & mpos_up);
    extern std::vector<Eigen::Tensor<cplx, 4>>  merge_unitary_mpo_layers(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> & mpos);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_time_evolution_operators_2site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &twosite_hams);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_time_evolution_operators_3site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &hams_3site);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_time_evolution_mpos(cplx delta_t, const std::vector<Eigen::Tensor<cplx, 4>> &mpos);
    extern cplx                                 get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &rho, size_t pos_rho, const Eigen::Matrix2cd &sig, size_t pos_sig);
    extern cplx                                 get_lbit_exp_value2(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern cplx                                 get_lbit_exp_value3(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern cplx                                 get_lbit_exp_value4(const std::vector<Eigen::Tensor<cplx, 4>> &mpo_layer, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj);
//    extern cplx                                 get_lbit_exp_value4(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern Eigen::Tensor<cplx, 2>               get_lbit_support(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites);
    extern Eigen::Tensor<cplx, 2>               get_lbit_support(const std::vector<Eigen::Tensor<cplx, 4>> &mpo_layer);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_lbit_supports(const UnitaryGateProperties &uprop, size_t reps, bool randomize_fields);
    extern std::pair<Eigen::Tensor<double, 2>,Eigen::Tensor<double, 2>>
                                                get_lbit_support_stats(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_support_vec);
    extern std::pair<Eigen::Tensor<double, 2>,Eigen::Tensor<double, 2>>
                                                get_lbit_permute_stats(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_support_vec);
    extern std::tuple<double,
                      double,
                      std::vector<double>,
                      size_t>                   get_characteristic_length_scale(const Eigen::Tensor<double, 2> &lbit_overlap_permuted);
//    extern lbitSupportAnalysis                  get_lbit_support_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec,
//                                                                          size_t reps, const UnitaryGateProperties & uprop, bool randomize_fields = false);
//    extern lbitSupportAnalysis                  get_lbit_support_analysis(const std::vector<UnitaryGateProperties> &uprops, bool randomize_fields = false);

    extern lbitSupportAnalysis                  get_lbit_support_analysis(
                                                                  const UnitaryGateProperties      u_defaults,
                                                                  size_t                           reps = 1,
                                                                  bool                             randomize_fields = false,
                                                                  std::vector<size_t          >    u_depths = {},
                                                                  std::vector<double          >    u_fmixs = {},
                                                                  std::vector<double          >    u_tstds = {},
                                                                  std::vector<double          >    u_cstds = {},
                                                                  std::vector<UnitaryGateWeight >  u_tgw8s = {},
                                                                  std::vector<UnitaryGateWeight >  u_cgw8s = {});

    /* clang-format on */
}