#pragma once
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
class StateInfinite;
class StateFinite;
class class_mps_2site;

namespace tools::common::views {
    using Scalar = std::complex<double>;
    extern Eigen::Tensor<Scalar, 4> theta, theta_evn_normalized, theta_odd_normalized;
    extern Eigen::Tensor<Scalar, 4> theta_sw;
    extern Eigen::Tensor<Scalar, 3> LAGA, LCGB;
    extern Eigen::Tensor<Scalar, 2> l_evn, r_evn;
    extern Eigen::Tensor<Scalar, 2> l_odd, r_odd;
    extern Eigen::Tensor<Scalar, 4> transfer_matrix_LAGA;
    extern Eigen::Tensor<Scalar, 4> transfer_matrix_LCGB;
    extern Eigen::Tensor<Scalar, 4> transfer_matrix_evn;
    extern Eigen::Tensor<Scalar, 4> transfer_matrix_odd;
    extern bool                     components_computed;
    extern void                     compute_mps_components(const StateInfinite &state);

    extern Eigen::Tensor<Scalar, 4> get_theta(const StateFinite &state, Scalar norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$.*/

    extern Eigen::Tensor<Scalar, 4> get_theta(const StateInfinite &state, Scalar norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    extern Eigen::Tensor<Scalar, 4> get_theta_swapped(const StateInfinite &state,
                                                      Scalar               norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    extern Eigen::Tensor<Scalar, 4> get_theta_evn(const StateInfinite &state, Scalar norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    extern Eigen::Tensor<Scalar, 4> get_theta_odd(const StateInfinite &state,
                                                  Scalar               norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_zero(const StateInfinite &state);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_LBGA(const StateInfinite &state, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_GALC(const StateInfinite &state, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_GBLB(const StateInfinite &state, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_LCGB(const StateInfinite &state, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_theta_evn(const StateInfinite &state, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_theta_odd(const StateInfinite &state, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_AB(const StateInfinite &state, int p);

    extern Eigen::Tensor<Scalar, 4> get_theta(const class_mps_2site &MPS, Scalar norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    extern Eigen::Tensor<Scalar, 4> get_theta_swapped(const class_mps_2site &MPS,
                                                      Scalar                 norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    extern Eigen::Tensor<Scalar, 4> get_theta_evn(const class_mps_2site &MPS, Scalar norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    extern Eigen::Tensor<Scalar, 4> get_theta_odd(const class_mps_2site &MPS,
                                                  Scalar                 norm = 1.0); /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_zero(const class_mps_2site &MPS);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_LBGA(const class_mps_2site &MPS, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_GALC(const class_mps_2site &MPS, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_GBLB(const class_mps_2site &MPS, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_LCGB(const class_mps_2site &MPS, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_theta_evn(const class_mps_2site &MPS, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_theta_odd(const class_mps_2site &MPS, Scalar norm = 1.0);
    extern Eigen::Tensor<Scalar, 4> get_transfer_matrix_AB(const class_mps_2site &MPS, int p);

}