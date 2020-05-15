#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <tools/common/log.h>
#include <tools/common/moments.h>
#include <config/nmspc_settings.h>

using cplx = std::complex<double>;
using real = double;
template<typename Scalar> using T3 = Eigen::Tensor<Scalar,3>;
template<typename Scalar> using T4 = Eigen::Tensor<Scalar,4>;
template <typename T> using TM = Eigen::TensorMap<T>;

template<typename mps_type, typename mpo_type, typename env_type>
extern double tools::common::moments::first(const mps_type &mps,const mpo_type &mpo, const env_type &envL, const env_type &envR){
    // This measures the expectation value <M> of some multisite operator M given multisite mps', mpos and corresponding environments.
    // This is usually the energy E = <H>
    // Note that the environments must contain the correct type of mpos
    using Scalar = typename mps_type::Scalar;
    static_assert(mpo_type::NumIndices == 4 and "Wrong mps tensor rank != 4 passed to calculation of first moment");
    static_assert(mpo_type::NumIndices == 4 and "Wrong mpo tensor rank != 4 passed to calculation of first moment");
    static_assert(env_type::NumIndices == 3 and "Wrong env tensor rank != 3 passed to calculation of first moment");

    Eigen::Tensor<Scalar, 0> M1;
    OMP omp(settings::threading::num_threads);
    /* clang-format off */
    M1.device(omp.dev) =
            envL
            .contract(mps,                  Textra::idx({0}, {1}))
            .contract(mpo,                  Textra::idx({2, 1}, {2, 0}))
            .contract(mps.conjugate(),      Textra::idx({3, 0}, {0, 1}))
            .contract(envR,                 Textra::idx({0, 2, 1}, {0, 1, 2}));
    /* clang-format on */

    if(abs(std::imag(M1(0))) > 1e-10) {
        tools::log->critical(fmt::format("First moment has an imaginary part: {:.16f} + i {:.16f}", std::real(M1(0)), std::imag(M1(0))));
        throw std::runtime_error(fmt::format("First moment has an imaginary part: {:.16f} + i {:.16f}", std::real(M1(0)), std::imag(M1(0))));
    }
    double moment = std::real(M1(0));
    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("First moment is invalid: {}", moment));
    return moment;
}


template double tools::common::moments::first(const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
template double tools::common::moments::first(const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
template double tools::common::moments::first(const TM<T3<real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
template double tools::common::moments::first(const TM<T3<cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);


template<typename mps_type, typename mpo_type, typename env_type>
extern double tools::common::moments::second(const mps_type &mps,const mpo_type &mpo, const env_type &envL, const env_type &envR) {
    // This measures the second moment <M²> of some multisite operator M given multisite mps', mpos and corresponding environments.
    // This is usually the second moment of the hamiltonian <H²>
    // Note that the environments must contain the correct type of mpos
    using Scalar = typename mps_type::Scalar;
    double log2chiL = std::log2(mps.dimension(1));
    double log2chiR = std::log2(mps.dimension(2));
    double log2spin = std::log2(mps.dimension(0));
    Eigen::Tensor<Scalar, 0> M2;
    OMP                      omp(settings::threading::num_threads);
    /* clang-format off */
    if(log2spin >= std::max(log2chiL, log2chiR)) {
        if(log2chiL > log2chiR) {
            //            tools::log->trace("H2 path: log2spin > std::max(log2chiL , log2chiR)  and  log2chiL > log2chiR ");
            Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
            M2.device(omp.dev) =
                mps_shuffled
                    .contract(envL, Textra::idx({0}, {0}))
                    .contract(mpo,  Textra::idx({0, 3}, {2, 0}))
                    .contract(envR, Textra::idx({0, 3}, {0, 2}))
                    .contract(mpo , Textra::idx({2, 1, 4}, {2, 0, 1}))
                    .contract(mps_shuffled.conjugate(), Textra::idx({2, 0, 1}, {1, 0, 2}));
        }
        else {
            //            tools::log->trace("H2 path: log2spin >= std::max(log2chiL , log2chiR) and  log2chiL <= log2chiR ");
            Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{2, 0, 1});
            M2.device(omp.dev) =
                mps_shuffled
                    .contract(envR, Textra::idx({0}, {0}))
                    .contract(mpo,  Textra::idx({0, 3}, {2, 1}))
                    .contract(envL, Textra::idx({0, 3}, {0, 2}))
                    .contract(mpo,  Textra::idx({2, 4, 1}, {2, 0, 1}))
                    .contract(mps_shuffled.conjugate(), Textra::idx({2, 1, 0}, {1, 2, 0}));
        }
    } else {
        //        tools::log->trace("H2 path: log2spin < std::max(log2chiL , log2chiR)");
        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
        M2.device(omp.dev) =
            mps_shuffled.contract(envL, Textra::idx({0}, {0}))
                .contract(mpo,  Textra::idx({0, 3}, {2, 0}))
                .contract(mpo,  Textra::idx({4, 2}, {2, 0}))
                .contract(envR, Textra::idx({0, 2, 3}, {0, 2, 3}))
                .contract(mps_shuffled.conjugate(), Textra::idx({1, 0, 2}, {1, 0, 2}));
    }
    /* clang-format on */
    if(abs(std::imag(M2(0))) > 1e-10) {
        tools::log->critical(fmt::format("Second moment has an imaginary part: {:.16f} + i {:.16f}", std::real(M2(0)), std::imag(M2(0))));
        throw std::runtime_error(fmt::format("Second moment has an imaginary part: {:.16f} + i {:.16f}", std::real(M2(0)), std::imag(M2(0))));
    }
    double moment = std::real(M2(0));
    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("Second moment is invalid: {}", moment));
    return moment;
}

template double tools::common::moments::second(const T3<real> &, const T4<real> &, const T4<real> &, const T4<real> &);
template double tools::common::moments::second(const T3<cplx> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);
template double tools::common::moments::second(const TM<T3<real>> &, const T4<real> &, const T4<real> &, const T4<real> &);
template double tools::common::moments::second(const TM<T3<cplx>> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);