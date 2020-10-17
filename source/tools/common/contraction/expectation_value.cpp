#include <general/nmspc_tensor_extra.h>
#include <tools/common/contraction.h>
#include <tools/common/fmt.h>
using namespace tools::common::contraction;

/* clang-format off */
template<typename Scalar>
double tools::common::contraction::expectation_value(const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                           const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                           const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                           const Scalar * const envR_ptr, std::array<long,3> envR_dims){

    // This measures the expectation value of some multisite mps with respect to some mpo operator and corresponding environments.
    // This is usually the energy E = <psi|H|psi> or variance V = <psi|(H-E)²|psi>
    // Note that the environments must contain the correct type of mpos
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps_ptr,mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(mpo_ptr,mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envR_ptr,envR_dims);

    if(mps.dimension(1) != envL.dimension(0))
        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions()));
    if(mps.dimension(2) != envR.dimension(0))
        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions()));
    if(mps.dimension(0) != mpo.dimension(2)) throw std::runtime_error(fmt::format("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions()));
    if(envL.dimension(2) != mpo.dimension(0))
        throw std::runtime_error(fmt::format("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions()));
    if(envR.dimension(2) != mpo.dimension(1))
        throw std::runtime_error(fmt::format("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions()));

    Eigen::Tensor<Scalar, 0> expval;
    expval.device(*Textra::omp::dev) =
        envL
            .contract(mps,                  Textra::idx({0}, {1}))
            .contract(mpo,                  Textra::idx({2, 1}, {2, 0}))
            .contract(mps.conjugate(),      Textra::idx({3, 0}, {0, 1}))
            .contract(envR,                 Textra::idx({0, 2, 1}, {0, 1, 2}));

    double moment = 0;
    if constexpr(std::is_same_v<Scalar,cplx>){
        moment = std::real(expval(0));
        if(abs(std::imag(expval(0))) > 1e-10)
            fmt::print("Expectation value has an imaginary part: {:.16f} + i {:.16f}\n", std::real(expval(0)), std::imag(expval(0)));
        //        throw std::runtime_error(fmt::format("Expectation value has an imaginary part: {:.16f} + i {:.16f}", std::real(expval(0)), std::imag(expval(0))));
    }else
        moment = expval(0);

    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("First moment is invalid: {}", moment));
    return moment;
}

template double tools::common::contraction::expectation_value(const real * const mps_ptr,  std::array<long,3> mps_dims,
                                                              const real * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                              const real * const envL_ptr, std::array<long,3> envL_dims,
                                                              const real * const envR_ptr, std::array<long,3> envR_dims);
template double tools::common::contraction::expectation_value(const cplx * const mps_ptr,  std::array<long,3> mps_dims,
                                                              const cplx * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                              const cplx * const envL_ptr, std::array<long,3> envL_dims,
                                                              const cplx * const envR_ptr, std::array<long,3> envR_dims);
/* clang-format on */



//
//template<typename mps_type, typename mpo_type, typename env_type>
//extern double tools::common::contraction::expectation_value(const mps_type &mps, const mpo_type &mpo, const env_type &envL, const env_type &envR) {
//    // This measures the expectation value <M> of some multisite operator M given multisite mps', mpos and corresponding environments.
//    // This is usually the energy E = <H>
//    // Note that the environments must contain the correct type of mpos
//    using Scalar = typename mps_type::Scalar;
//    static_assert(mpo_type::NumIndices == 4 and "Wrong mps tensor rank != 4 passed to calculation of expectation_value");
//    static_assert(mpo_type::NumIndices == 4 and "Wrong mpo tensor rank != 4 passed to calculation of expectation_value");
//    static_assert(env_type::NumIndices == 3 and "Wrong env tensor rank != 3 passed to calculation of expectation_value");
//
//    if(mps.dimension(1) != envL.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions()));
//    if(mps.dimension(2) != envR.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions()));
//    if(mps.dimension(0) != mpo.dimension(2)) throw std::runtime_error(fmt::format("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions()));
//
//    Eigen::Tensor<Scalar, 0> M1;
//    /* clang-format off */
//    M1.device(*Textra::omp::dev) =
//        envL
//            .contract(mps,                  Textra::idx({0}, {1}))
//            .contract(mpo,                  Textra::idx({2, 1}, {2, 0}))
//            .contract(mps.conjugate(),      Textra::idx({3, 0}, {0, 1}))
//            .contract(envR,                 Textra::idx({0, 2, 1}, {0, 1, 2}));
//    /* clang-format on */
//
//    if(abs(std::imag(M1(0))) > 1e-10) {
//        throw std::runtime_error(fmt::format("First moment has an imaginary part: {:.16f} + i {:.16f}", std::real(M1(0)), std::imag(M1(0))));
//    }
//    double moment = std::real(M1(0));
//    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("First moment is invalid: {}", moment));
//    return moment;
//}
//
//
//template double tools::common::contraction::expectation_value(const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
//template double tools::common::contraction::expectation_value(const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
//template double tools::common::contraction::expectation_value(const TM<T3<real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
//template double tools::common::contraction::expectation_value(const TM<T3<cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);

