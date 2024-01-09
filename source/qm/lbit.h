#pragma once
#include "gate.h"
#include "math/svd/config.h"
#include "math/tenx/fwd_decl.h"
#include "qm.h"
#include <complex>
#include <optional>
#include <string_view>
#include <vector>

enum class UnitaryGateWeight;
enum class UnitaryGateType;
enum class MeanType;
class StateFinite;
namespace h5pp {
    class File;
    namespace hid {
        class h5t;
    }
};

namespace qm::lbit {

    /*! The values that we need to recreate a unitary gate exactly */
    struct UnitaryGateParameters {
        size_t                layer; /*!< Layer number of this gate */
        std::array<size_t, 2> sites; /*!< The two sites which this site connects */
        double                f;     /*!< Mixing factor  */
        double                w;     /*!< Gate weight value = exp(-2|h_i - h_i+1|) */
        std::array<double, 4> theta; /*!< Random real valued thetas (for the mixing term)  */
        std::complex<double>  c;     /*!< Random complex valued c (exchange term)  */
        UnitaryGateType       type;  /*!< Gate type non-interacting or interacting: ANDERSON(0) OR MBL(1)  */

        [[nodiscard]] static const h5pp::hid::h5t          get_h5_type();
        [[nodiscard]] static std::vector<std::string_view> get_parameter_names() noexcept;
        void                                               print_parameter_names() noexcept;
        void                                               print_parameter_values() const noexcept;
        [[nodiscard]] std::string                          fmt_value(std::string_view p) const;

        private:
        static h5pp::hid::h5t get_h5t_enum_ut();
    };

    struct UnitaryGateProperties {
        size_t                                     sites;   /*!< Width of the circuit, i.e. system length/number of sites */
        size_t                                     depth;   /*!< Number of layers in the unitary circuit (1 layer connects all neighbors once) */
        double                                     fmix;    /*!< The mixing factor f in exp(-ifM) */
        double                                     tstd;    /*!< Standard deviation for the theta parameters in the unitary gate */
        double                                     cstd;    /*!< Standard deviation for the c parameters in the unitary gate  */
        UnitaryGateWeight                          g8w8;    /*!< Choose IDENTITY|EXPDECAY type of gate weight (hvals is required for EXPDECAY) */
        UnitaryGateType                            type;    /*!< Choose ANDERSON|MBL */
        double                                     hmean;   /*!< mean of random onsite fields */
        double                                     hwdth;   /*!< width of random onsite fields (st.dev. if normal) */
        std::string_view                           hdist;   /*!< distribution of onsite fields */
        mutable std::vector<double>                hvals;   /*!< onsite fields of the l-bit hamiltonian, needed for type == UnitaryGateWeight::EXPDECAY */
        mutable std::vector<std::vector<qm::Gate>> ulayers; /*!< The generated unitary circuit of two-site gates */
        UnitaryGateProperties() = default;
        UnitaryGateProperties(const std::vector<double> &h);
        std::string    string() const;
        void           randomize_hvals() const;
        bool           keep_circuit = false;
        mutable size_t layer_count  = 0;

        mutable std::vector<UnitaryGateParameters> circuit;
    };

    void write_unitary_circuit_parameters(h5pp::File &file, std::string_view table_path, const std::vector<UnitaryGateParameters> &circuit);
    std::vector<UnitaryGateParameters> read_unitary_circuit_parameters(const h5pp::File &file, std::string_view table_path);
    std::vector<std::vector<qm::Gate>> read_unitary_2site_gate_layers(const h5pp::File &file, std::string_view table_path);
    std::vector<std::vector<qm::Gate>> get_unitary_2site_gate_layers(const std::vector<UnitaryGateParameters> &circuit);

    struct lbitSupportAnalysis {
        Eigen::Tensor<real, 6> cls_avg_fit; // Characteristic length-scale of lbits from linear regression of log data
        Eigen::Tensor<real, 6> cls_avg_rms; // Root mean squared deviation
        Eigen::Tensor<real, 6> cls_avg_rsq; // R-squared or coefficient of determination
        Eigen::Tensor<real, 6> cls_typ_fit; // Characteristic length-scale of lbits from linear regression of log data
        Eigen::Tensor<real, 6> cls_typ_rms; // Root mean squared deviation
        Eigen::Tensor<real, 6> cls_typ_rsq; // R-squared or coefficient of determination
        Eigen::Tensor<real, 7> corravg;     // The lbit correlation matrix of l-bits averaged over site and disorder
        Eigen::Tensor<real, 7> corrtyp;     // The lbit correlation matrix of l-bits geometrically averaged over site and disorder
        Eigen::Tensor<real, 7> correrr;     // The sterr of l-bits: permuted and averaged over site and disorder
        Eigen::Tensor<real, 9> corrmat;     // The raw data from l-bit correlation matrices or the trace O(i,j) for each realization
        Eigen::Tensor<real, 9> corroff;     // The offset lbit correlation matrices O(i, |i-j|) for each realization.
        lbitSupportAnalysis() {
            cls_avg_fit.setZero();
            cls_avg_rms.setZero();
            cls_avg_rsq.setZero();
            cls_typ_fit.setZero();
            cls_typ_rms.setZero();
            cls_typ_rsq.setZero();
            corravg.setZero();
            corrtyp.setZero();
            correrr.setZero();
            corrmat.setZero();
            corroff.setZero();
        }
        lbitSupportAnalysis(size_t ndpth, size_t nfmix, size_t ntstd, size_t ncstd, size_t ng8w8, size_t ntype, size_t nreps, size_t nsize)
            : lbitSupportAnalysis() {
            auto idpth = static_cast<Eigen::Index>(ndpth);
            auto ifmix = static_cast<Eigen::Index>(nfmix);
            auto itstd = static_cast<Eigen::Index>(ntstd);
            auto icstd = static_cast<Eigen::Index>(ncstd);
            auto ig8w8 = static_cast<Eigen::Index>(ng8w8);
            auto itype = static_cast<Eigen::Index>(ntype);
            auto ireps = static_cast<Eigen::Index>(nreps);
            auto isize = static_cast<Eigen::Index>(nsize);
            cls_avg_fit.resize(idpth, ifmix, itstd, icstd, ig8w8, itype);
            cls_avg_rms.resize(idpth, ifmix, itstd, icstd, ig8w8, itype);
            cls_avg_rsq.resize(idpth, ifmix, itstd, icstd, ig8w8, itype);
            cls_typ_fit.resize(idpth, ifmix, itstd, icstd, ig8w8, itype);
            cls_typ_rms.resize(idpth, ifmix, itstd, icstd, ig8w8, itype);
            cls_typ_rsq.resize(idpth, ifmix, itstd, icstd, ig8w8, itype);
            corravg.resize(idpth, ifmix, itstd, icstd, ig8w8, itype, isize);
            corrtyp.resize(idpth, ifmix, itstd, icstd, ig8w8, itype, isize);
            correrr.resize(idpth, ifmix, itstd, icstd, ig8w8, itype, isize);
            corrmat.resize(idpth, ifmix, itstd, icstd, ig8w8, itype, ireps, isize, isize);
            corroff.resize(idpth, ifmix, itstd, icstd, ig8w8, itype, ireps, isize, isize);
        }
    };

    /* clang-format off */
    extern Eigen::Tensor<cplx, 2>               get_unitary_layer_as_tensor(const std::vector<qm::Gate> &unitary_layer);
    extern Eigen::Tensor<cplx, 2>               get_unitary_circuit_as_tensor(const std::vector<std::vector<qm::Gate>> &unitary_circuit);
    extern Eigen::Tensor<cplx, 2>               get_time_evolution_operator(cplx_t delta_t, const Eigen::Tensor<cplx, 2> &hamiltonian);
    extern std::vector<qm::Gate>                get_time_evolution_gates(cplx_t delta_t, const std::vector<qm::Gate> &hams_nsite, double id_threshold = std::numeric_limits<double>::epsilon());
    extern std::vector<qm::SwapGate>            get_time_evolution_swap_gates(cplx_t delta_t, const std::vector<qm::SwapGate> &hams_nsite, double id_threshold = std::numeric_limits<double>::epsilon());
//    extern std::vector<qm::Gate>                get_unitary_2gate_layer(size_t sites, double fmix);
    extern qm::Gate                             get_unitary_2site_gate(const UnitaryGateParameters &u);
    extern std::vector<qm::Gate>                create_unitary_2site_gate_layer(const qm::lbit::UnitaryGateProperties &u);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_unitary_mpo_layer(const std::vector<qm::Gate> & ulayer, std::optional<svd::config> cfg = std::nullopt);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_unitary_mpo_layer(const UnitaryGateProperties & u);
    extern std::vector<Eigen::Tensor<cplx, 4>>  merge_unitary_mpo_layers(const std::vector<Eigen::Tensor<cplx, 4>> & mpos_dn, const std::vector<Eigen::Tensor<cplx, 4>> & mpos_up, bool adj_dn = false);
    extern std::vector<Eigen::Tensor<cplx, 4>>  merge_unitary_mpo_layers(const std::vector<Eigen::Tensor<cplx, 4>> & mpos_dn,
                                                                         const std::vector<Eigen::Tensor<cplx, 4>> & mpos_md,
                                                                         const std::vector<Eigen::Tensor<cplx, 4>> & mpos_up);
    extern std::vector<Eigen::Tensor<cplx, 4>>  merge_unitary_mpo_layers(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> & mpos);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_time_evolution_operators_2site(size_t sites, cplx_t delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &twosite_hams);
    extern std::vector<Eigen::Tensor<cplx, 2>>  get_time_evolution_operators_3site(size_t sites, cplx_t delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &hams_3site);
    extern std::vector<Eigen::Tensor<cplx, 4>>  get_time_evolution_mpos(cplx_t delta_t, const std::vector<Eigen::Tensor<cplx, 4>> &mpos);
    extern cplx                                 get_lbit_2point_correlator1(const std::vector<std::vector<qm::Gate>> &unitary_circuit, const Eigen::Matrix2cd &rho, size_t pos_rho, const Eigen::Matrix2cd &sig, size_t pos_sig);
    extern cplx                                 get_lbit_2point_correlator2(const std::vector<std::vector<qm::Gate>> &unitary_circuit, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern cplx                                 get_lbit_2point_correlator3(const std::vector<std::vector<qm::Gate>> &unitary_circuit, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj, long len);
    extern cplx                                 get_lbit_2point_correlator4(const std::vector<Eigen::Tensor<cplx, 4>> &mpo_layer, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj, size_t pos_szj);
    extern Eigen::Tensor<cplx, 1>               get_lbit_2point_correlator5(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &mpo_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj);
    extern Eigen::Tensor<cplx, 1>               get_lbit_2point_correlator6(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &mpo_layers, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj);
    extern Eigen::Tensor<real, 2>               get_lbit_correlation_matrix(const std::vector<Eigen::Tensor<cplx, 4>> &mpo_layer);
    extern Eigen::Tensor<real, 2>               get_lbit_correlation_matrix(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &mpo_layers);
    extern Eigen::Tensor<real, 2>               get_lbit_correlation_matrix2(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &mpo_layers);
    extern Eigen::Tensor<real, 2>               get_lbit_correlation_matrix(const std::vector<std::vector<qm::Gate>> &unitary_circuit, size_t sites);
    extern Eigen::Tensor<real, 2>               get_lbit_correlation_matrix(const std::vector<std::vector<qm::Gate>> &unitary_circuit, size_t sites, size_t max_num_states, double tol);
    extern std::vector<Eigen::Tensor<real, 2>>  get_lbit_correlation_matrices(const UnitaryGateProperties &uprop, size_t reps, bool randomize_fields, bool use_mpo);
    extern std::vector<Eigen::Tensor<real, 2>>  get_lbit_correlation_matrices2(const UnitaryGateProperties &uprop, size_t reps, bool randomize_fields, bool use_mpo);
    extern std::tuple<Eigen::Tensor<real, 2>, Eigen::Tensor<real, 2>, Eigen::Tensor<real, 2>>
                                                get_lbit_correlation_statistics(const std::vector<Eigen::Tensor<real, 2>> &lbit_corrmats);
    extern std::tuple<double,
                      double,
                      double,
                      std::vector<double>,
                      size_t>                   get_characteristic_length_scale(const Eigen::Tensor<real, 2> &lbit_corrmat_disorder_mean, MeanType);
//    extern lbitSupportAnalysis                  get_lbit_support_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec,
//                                                                          size_t reps, const UnitaryGateProperties & uprop, bool randomize_fields = false);
//    extern lbitSupportAnalysis                  get_lbit_support_analysis(const std::vector<UnitaryGateProperties> &uprops, bool randomize_fields = false);

    extern lbitSupportAnalysis                  get_lbit_support_analysis(
                                                                  const UnitaryGateProperties      & u_defaults,
                                                                  std::vector<size_t          >    udpths = {},
                                                                  std::vector<double          >    ufmixs = {},
                                                                  std::vector<double          >    utstds = {},
                                                                  std::vector<double          >    ucstds = {},
                                                                  std::vector<UnitaryGateWeight >  ug8w8s = {},
                                                                  std::vector<UnitaryGateType   >  utypes = {}
                                                                );

    StateFinite transform_to_real_basis(const StateFinite &lbit_state,
                                        const std::vector<std::vector<qm::Gate>> & unitary_gates_2site_layers,
                                        svd::config svd_cfg);
    StateFinite transform_to_real_basis(const StateFinite &lbit_state,
                                        const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &unitary_gates_mpo_layers,
                                        const Eigen::Tensor<std::complex<double>, 1> & ledge,
                                        const Eigen::Tensor<std::complex<double>, 1> & redge,
                                        svd::config svd_cfg);
    StateFinite transform_to_lbit_basis(const StateFinite &real_state,
                                        const std::vector<std::vector<qm::Gate>> & unitary_gates_2site_layers,
                                        svd::config svd_cfg);
    StateFinite transform_to_lbit_basis(const StateFinite &real_state,
                                        const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &unitary_gates_mpo_layers,
                                        const Eigen::Tensor<std::complex<double>, 1> & ledge,
                                        const Eigen::Tensor<std::complex<double>, 1> & redge,
                                        svd::config svd_cfg);
    /* clang-format on */
}