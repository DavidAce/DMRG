#pragma once
#include <complex>
#include <deque>
#include <general/eigen_tensor_fwd_decl.h>
#include <math/svd/settings.h>
#include <optional>
#include <tuple>
#include <vector>
class MpsSite;
namespace tools::common::split {
    using real = double;
    using cplx = std::complex<double>;
    /* clang-format off */
    template<typename Scalar>
    extern std::vector<MpsSite> split_mps (const Eigen::Tensor<Scalar,3> & multisite_mps,
                                           const std::vector<long>         & spin_dims,
                                           const std::vector<size_t>       & positions,
                                           long                            center_position,
                                           long                            chi_limit,
                                           std::optional<svd::settings>    svd_settings = std::nullopt);


    namespace internal{
        template<typename Scalar>
        extern std::vector<MpsSite>
                    split_mps_into_As(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                      const std::vector<long>       & spin_dims,
                                      const std::vector<size_t>     & positions,
                                      long                            chi_limit,
                                      std::optional<svd::settings>    svd_settings = std::nullopt);
        template<typename Scalar>
        extern std::deque<MpsSite>
                    split_mps_into_Bs(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                      const std::vector<long>       & spin_dims,
                                      const std::vector<size_t>     & positions,
                                      long                            chi_limit,
                                      std::optional<svd::settings>    svd_settings = std::nullopt);
    }

    /* clang-format on */
}