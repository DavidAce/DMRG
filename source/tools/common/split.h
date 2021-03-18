#pragma once
#include <complex>
#include <deque>
#include <general/eigen_tensor_fwd_decl.h>
#include <optional>
#include <tuple>
#include <vector>
#include <math/svd/settings.h>
class class_mps_site;
namespace tools::common::split {
    using Scalar = std::complex<double>;
    /* clang-format off */
    extern std::vector<class_mps_site> split_mps (const Eigen::Tensor<Scalar,3> & multisite_mps,
                                                  const std::vector<long>         & spin_dims,
                                                  const std::vector<size_t>       & positions,
                                                  long                            center_position,
                                                  long                            chi_limit,
                                                  std::optional<svd::settings>    svd_settings = std::nullopt);


    namespace internal{

        extern std::vector<class_mps_site>
                    split_mps_into_As(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                      const std::vector<long>       & spin_dims,
                                      const std::vector<size_t>     & positions,
                                      long                            chi_limit,
                                      std::optional<svd::settings>    svd_settings = std::nullopt);

        extern std::deque<class_mps_site>
                    split_mps_into_Bs(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                      const std::vector<long>       & spin_dims,
                                      const std::vector<size_t>     & positions,
                                      long                            chi_limit,
                                      std::optional<svd::settings>    svd_settings = std::nullopt);
    }

    /* clang-format on */
}