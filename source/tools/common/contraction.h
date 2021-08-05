#pragma once

#include <general/eigen_tensor_fwd_decl.h>
// Eigen goes first

#include <complex>

namespace tools::common::contraction {
    using cplx = std::complex<double>;
    using real = double;

    /* clang-format off */
    template<typename Scalar>
    double expectation_value(  const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                               const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                               const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                               const Scalar * const envR_ptr, std::array<long,3> envR_dims);

    template<typename Scalar>
    void matrix_vector_product(Scalar * result,
                               const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                               const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                               const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                               const Scalar * const envR_ptr, std::array<long,3> envR_dims);

    template<typename Scalar>
    void matrix_inverse_vector_product(Scalar * result,
                               const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                               const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                               const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                               const Scalar * const envR_ptr, std::array<long,3> envR_dims);


    template<typename mps_type, typename mpo_type, typename env_type>
    double expectation_value(
        const Eigen::TensorBase<mps_type,Eigen::ReadOnlyAccessors> & mps,
        const Eigen::TensorBase<mpo_type,Eigen::ReadOnlyAccessors> & mpo,
        const Eigen::TensorBase<env_type,Eigen::ReadOnlyAccessors> & envL,
        const Eigen::TensorBase<env_type,Eigen::ReadOnlyAccessors> & envR){
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3 passed to calculation of expectation_value");
        static_assert(mpo_type::NumIndices == 4 and "Wrong mpo tensor rank != 4 passed to calculation of expectation_value");
        static_assert(env_type::NumIndices == 3 and "Wrong env tensor rank != 3 passed to calculation of expectation_value");
        const auto & mps_ref  = static_cast<const mps_type&>(mps);
        const auto & mpo_ref  = static_cast<const mpo_type&>(mpo);
        const auto & envL_ref = static_cast<const env_type&>(envL);
        const auto & envR_ref = static_cast<const env_type&>(envR);
        return expectation_value(
            mps_ref.data(), mps_ref.dimensions(),
            mpo_ref.data(), mpo_ref.dimensions(),
            envL_ref.data(), envL_ref.dimensions(),
            envR_ref.data(), envR_ref.dimensions());
    }

    template<typename res_type, typename mps_type, typename mpo_type, typename env_type>
    void matrix_vector_product(
        Eigen::TensorBase<res_type,Eigen::WriteAccessors> &res,
        const Eigen::TensorBase<mps_type,Eigen::ReadOnlyAccessors> & mps,
        const Eigen::TensorBase<mpo_type,Eigen::ReadOnlyAccessors> & mpo,
        const Eigen::TensorBase<env_type,Eigen::ReadOnlyAccessors> & envL,
        const Eigen::TensorBase<env_type,Eigen::ReadOnlyAccessors> & envR){
        static_assert(res_type::NumIndices == 3 and "Wrong res tensor rank != 3 passed to calculation of matrix_vector_product");
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3 passed to calculation of matrix_vector_product");
        static_assert(mpo_type::NumIndices == 4 and "Wrong mpo tensor rank != 4 passed to calculation of matrix_vector_product");
        static_assert(env_type::NumIndices == 3 and "Wrong env tensor rank != 3 passed to calculation of matrix_vector_product");
        auto & res_ref = static_cast<res_type&>(res);
        const auto & mps_ref  = static_cast<const mps_type&>(mps);
        const auto & mpo_ref  = static_cast<const mpo_type&>(mpo);
        const auto & envL_ref = static_cast<const env_type&>(envL);
        const auto & envR_ref = static_cast<const env_type&>(envR);
        matrix_vector_product(
            res_ref.data(),
            mps_ref.data(), mps_ref.dimensions(),
            mpo_ref.data(), mpo_ref.dimensions(),
            envL_ref.data(), envL_ref.dimensions(),
            envR_ref.data(), envR_ref.dimensions());
    }

    template<typename mps_type, typename mpo_type, typename env_type>
    mps_type matrix_vector_product(const mps_type &mps, const mpo_type &mpo, const env_type &envL, const env_type &envR) {
        mps_type result;
        result.resize(mps.dimensions());
        matrix_vector_product(result, mps, mpo, envL, envR);
        return result;
    }

    template<typename res_type, typename mps_type, typename mpo_type, typename env_type>
    void matrix_inverse_vector_product(
        Eigen::TensorBase<res_type,Eigen::WriteAccessors> &res,
        const Eigen::TensorBase<mps_type,Eigen::ReadOnlyAccessors> & mps,
        const Eigen::TensorBase<mpo_type,Eigen::ReadOnlyAccessors> & mpo,
        const Eigen::TensorBase<env_type,Eigen::ReadOnlyAccessors> & envL,
        const Eigen::TensorBase<env_type,Eigen::ReadOnlyAccessors> & envR){
        static_assert(res_type::NumIndices == 3 and "Wrong res tensor rank != 3 passed to calculation of matrix_vector_product");
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3 passed to calculation of matrix_vector_product");
        static_assert(mpo_type::NumIndices == 4 and "Wrong mpo tensor rank != 4 passed to calculation of matrix_vector_product");
        static_assert(env_type::NumIndices == 3 and "Wrong env tensor rank != 3 passed to calculation of matrix_vector_product");
        auto & res_ref = static_cast<res_type&>(res);
        const auto & mps_ref  = static_cast<const mps_type&>(mps);
        const auto & mpo_ref  = static_cast<const mpo_type&>(mpo);
        const auto & envL_ref = static_cast<const env_type&>(envL);
        const auto & envR_ref = static_cast<const env_type&>(envR);
        matrix_inverse_vector_product(
            res_ref.data(),
            mps_ref.data(), mps_ref.dimensions(),
            mpo_ref.data(), mpo_ref.dimensions(),
            envL_ref.data(), envL_ref.dimensions(),
            envR_ref.data(), envR_ref.dimensions());
    }


    /* clang-format on */

    // Extern templates

    //    template<typename Scalar>
    //    using T3 = Eigen::Tensor<Scalar, 3>;
    //    template<typename Scalar>
    //    using T4 = Eigen::Tensor<Scalar, 4>;
    //    template<typename T>
    //    using TM = Eigen::TensorMap<T>;
    //    extern template double expectation_value(const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
    //    extern template double expectation_value(const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);

    //    extern template double expectation_value(const TM<T3<real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
    //    extern template double expectation_value(const TM<T3<cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
    //    extern template double second_moment(const T3<real> &, const T4<real> &, const T4<real> &, const T4<real> &);
    //    extern template double second_moment(const T3<cplx> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);
    //    extern template double second_moment(const TM<T3<real>> &, const T4<real> &, const T4<real> &, const T4<real> &);
    //    extern template double second_moment(const TM<T3<cplx>> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);

    //    extern template T3<real> matrix_vector_product(const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
    //    extern template T3<cplx> matrix_vector_product(const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
    //    extern template void matrix_vector_product(T3<real> & ,const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
    //    extern template void matrix_vector_product(T3<cplx> & ,const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
    //    extern template void matrix_vector_product(T3<real> & ,const TM<T3<const real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
    //    extern template void matrix_vector_product(T3<cplx> & ,const TM<T3<const cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);

}