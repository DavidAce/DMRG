#pragma omp
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>

class class_state_finite;
class class_model_finite;

namespace tools::finite::views{
    extern const Eigen::Tensor<std::complex<double>,3> & get_multisite_mps(const class_state_finite & state, const std::list<size_t> & active_sites);
    extern const Eigen::Tensor<std::complex<double>,3> & get_multisite_mps(const class_state_finite & state);
    extern const Eigen::Tensor<std::complex<double>,4> & get_multisite_mpo(const class_model_finite & model, const std::list<size_t> & active_sites);
    extern const Eigen::Tensor<std::complex<double>,4> & get_multisite_mpo(const class_model_finite & model);

}