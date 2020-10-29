#pragma once
#include <complex>
#include <general/eigen_tensor_fwd_decl.h>
#include <list>
#include <tuple>
#include <vector>
#include <optional>
class class_mps_site;
namespace tools::common::split {
    using Scalar = std::complex<double>;
    /* clang-format off */
    extern std::list<class_mps_site> split_mps (const Eigen::Tensor<Scalar,3> & multisite_mps,
                                                const std::vector<long>         & spin_dims,
                                                const std::vector<size_t>       & sites,
                                                size_t                          center_position,
                                                long                            chi_limit,
                                                std::optional<double>           svd_threshold = std::nullopt);


    namespace internal{

        extern std::list<class_mps_site>
//            extern std::pair<std::list<class_mps_site>,Eigen::Tensor<Scalar,3>>
                    split_mps_from_left(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                        std::vector<long>               spin_dims,
                                        std::vector<size_t>             positions,
                                        long                            chi_limit,
                                        std::optional<double>           svd_threshold = std::nullopt);

        extern std::list<class_mps_site>
//            extern std::pair<Eigen::Tensor<Scalar,3>,std::list<class_mps_site>>
            split_mps_from_right(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                         std::vector<long>               spin_dims,
                                         std::vector<size_t>             positions,
                                         long                            chi_limit,
                                         std::optional<double>           svd_threshold = std::nullopt);
    }

    /* clang-format on */
}