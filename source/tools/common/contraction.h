#pragma once

#include <math/tenx/fwd_decl.h>
// Eigen goes first
#include <complex>
#include <debug/exceptions.h>
#include <math/tenx.h>

namespace tools::common::contraction {
    using cplx = std::complex<double>;
    using real = double;
    template<typename T>
    using TensorWrite = Eigen::TensorBase<T, Eigen::WriteAccessors>;
    template<typename T>
    using TensorRead = Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>;

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
        const TensorRead<mps_type> & mps,
        const TensorRead<mpo_type> & mpo,
        const TensorRead<env_type> & envL,
        const TensorRead<env_type> & envR){
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
        TensorWrite<res_type> &res,
        const TensorRead<mps_type> & mps,
        const TensorRead<mpo_type> & mpo,
        const TensorRead<env_type> & envL,
        const TensorRead<env_type> & envR){
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
        TensorWrite<res_type> &res,
        const TensorRead<mps_type> & mps,
        const TensorRead<mpo_type> & mpo,
        const TensorRead<env_type> & envL,
        const TensorRead<env_type> & envR){
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
    template<typename T>
    using OptRefWrap = std::optional<std::reference_wrapper<T>>;
    template<typename T>
    using TensorWriteOptRef = OptRefWrap<TensorWrite<T>>;

    template<typename Scalar>
    void contract_mps_bnd(Scalar *res_ptr, std::array<long, 3> res_dims, const Scalar *const mps_ptr, std::array<long, 3> mps_dims, const Scalar *const bnd_ptr,
                          std::array<long, 1> bnd_dims);

    template<typename Scalar>
    void contract_bnd_mps(Scalar *res_ptr, std::array<long, 3> res_dims, const Scalar *const bnd_ptr, std::array<long, 1> bnd_dims, const Scalar *const mps_ptr,
                          std::array<long, 3> mps_dims);

    template<typename res_type, typename mps_type, typename bnd_type>
    void contract_mps_bnd(TensorWrite<res_type> &res, const TensorRead<mps_type> &mps, const TensorRead<bnd_type> &bnd) {
        static_assert(res_type::NumIndices == 3 and "Wrong res tensor rank != 3");
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        static_assert(bnd_type::NumIndices == 1 and "Wrong bnd tensor rank != 1");
        auto &res_ref  = static_cast<res_type &>(res);
        auto  mps_eval = tenx::asEval(mps, tenx::omp::getDevice());
        auto  bnd_eval = tenx::asEval(bnd, tenx::omp::getDevice());

        res_ref.resize(mps_eval.dimensions());
        contract_mps_bnd(res_ref.data(), res_ref.dimensions(), mps_eval.data(), mps_eval.dimensions(), bnd_eval.data(), bnd_eval.dimensions());
    }

    template<typename res_type, typename bnd_type, typename mps_type>
    void contract_bnd_mps(TensorWrite<res_type> &res, const TensorRead<bnd_type> &bnd, const TensorRead<mps_type> &mps) {
        static_assert(res_type::NumIndices == 3 and "Wrong res tensor rank != 3");
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        static_assert(bnd_type::NumIndices == 1 and "Wrong bnd tensor rank != 1");
        auto       &res_ref = static_cast<res_type &>(res);
        const auto &mps_ref = static_cast<const mps_type &>(mps);
        const auto &bnd_ref = static_cast<const bnd_type &>(bnd);
        res_ref.resize(mps_ref.dimensions());
        contract_bnd_mps(res_ref.data(), res_ref.dimensions(), bnd_ref.data(), bnd_ref.dimensions(), mps_ref.data(), mps_ref.dimensions());
    }

    template<typename mps_type, typename bnd_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 3> contract_mps_bnd_temp(const TensorRead<mps_type> &mps, const TensorRead<bnd_type> &bnd,
                                                                                    Eigen::Tensor<typename mps_type::Scalar, 3> &temp) {
        contract_mps_bnd(temp, mps, bnd);
        return temp;
    }
    template<typename mps_type, typename bnd_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 3> contract_mps_bnd_temp(const TensorRead<mps_type> &mps, const TensorRead<bnd_type> &bnd) {
        Eigen::Tensor<typename mps_type::Scalar, 3> temp;
        return contract_mps_bnd_temp(mps, bnd, temp);
    }

    template<typename bnd_type, typename mps_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 3> contract_bnd_mps_temp(const TensorRead<bnd_type> &bnd, const TensorRead<mps_type> &mps,
                                                                                    Eigen::Tensor<typename mps_type::Scalar, 3> &temp) {
        contract_bnd_mps(temp, bnd, mps);
        return temp;
    }

    template<typename bnd_type, typename mps_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 3> contract_bnd_mps_temp(const TensorRead<bnd_type> &bnd, const TensorRead<mps_type> &mps) {
        Eigen::Tensor<typename mps_type::Scalar, 3> temp;
        return contract_bnd_mps_temp(bnd, mps, temp);
    }

    template<typename Scalar>
    void contract_mps_mps(Scalar *res_ptr, std::array<long, 3> res_dims, const Scalar *const mpsL_ptr, std::array<long, 3> mpsL_dims,
                          const Scalar *const mpsR_ptr, std::array<long, 3> mpsR_dims);

    template<typename mps_type>
    void contract_mps_mps(TensorWrite<mps_type> &res, const TensorRead<mps_type> &mpsL, const TensorRead<mps_type> &mpsR) {
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        auto               &res_ref  = static_cast<mps_type &>(res);
        const auto         &mpsL_ref = static_cast<const mps_type &>(mpsL);
        const auto         &mpsR_ref = static_cast<const mps_type &>(mpsR);
        long                d0       = mpsL_ref.dimension(0) * mpsR_ref.dimension(0);
        long                d1       = mpsL_ref.dimension(1);
        long                d2       = mpsR_ref.dimension(2);
        std::array<long, 3> res_dims = {d0, d1, d2};
        res_ref.resize(res_dims);
        contract_mps_mps(res_ref.data(), res_ref.dimensions(), mpsL_ref.data(), mpsL_ref.dimensions(), mpsR_ref.data(), mpsR_ref.dimensions());
    }

    template<typename mps_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 3> contract_mps_mps_temp(const TensorRead<mps_type> &mpsL, const TensorRead<mps_type> &mpsR,
                                                                                    Eigen::Tensor<typename mps_type::Scalar, 3> &temp) {
        contract_mps_mps(temp, mpsL, mpsR);
        return temp;
    }

    template<typename mps_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 3> contract_mps_mps_temp(const TensorRead<mps_type> &mpsL, const TensorRead<mps_type> &mpsR) {
        Eigen::Tensor<typename mps_type::Scalar, 3> temp;
        return contract_mps_mps_temp(mpsL, mpsR, temp);
    }

    template<typename Scalar>
    double contract_mps_mps_overlap(const Scalar *const mps1_ptr, std::array<long, 3> mps1_dims, const Scalar *const mps2_ptr, std::array<long, 3> mps2_dims);

    template<typename mps_type>
    double contract_mps_norm(const TensorRead<mps_type> &mps1, const TensorRead<mps_type> &mps2) {
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        const auto &mps1_ref = static_cast<const mps_type &>(mps1);
        const auto &mps2_ref = static_cast<const mps_type &>(mps2);
        return contract_mps_overlap(mps1_ref.data(), mps1_ref.dimensions(), mps2_ref.data(), mps2_ref.dimensions());
    }

    template<typename mps_type>
    double contract_mps_norm(const TensorRead<mps_type> &mps) {
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        const auto &mps_ref = static_cast<const mps_type &>(mps);
        return contract_mps_mps_overlap(mps_ref.data(), mps_ref.dimensions(), mps_ref.data(), mps_ref.dimensions());
    }

    template<typename Scalar>
    void contract_mps_mps_partial(Scalar *res_ptr, std::array<long, 2> res_dims, const Scalar *const mps1_ptr, std::array<long, 3> mps1_dims,
                                  const Scalar *const mps2_ptr, std::array<long, 3> mps2_dims, std::array<long, 2> idx);

    template<typename mps_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 2> contract_mps_mps_partial(const TensorRead<mps_type> &mps1, const TensorRead<mps_type> &mps2,
                                                                                       std::array<long, 2> idx) {
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        const auto         &mps1_ref = static_cast<const mps_type &>(mps1);
        const auto         &mps2_ref = static_cast<const mps_type &>(mps2);
        std::array<long, 2> dims     = {};
        if(idx[1] == 1)
            dims = {mps1_ref.dimension(2), mps2_ref.dimension(2)};
        else if(idx[1] == 2)
            dims = {mps1_ref.dimension(1), mps2_ref.dimension(1)};
        else
            throw except::runtime_error("Expected idx [0,1] or [0,2]. Got: {}", idx);
        Eigen::Tensor<typename mps_type::Scalar, 2> res(dims);
        contract_mps_mps_partial(res.data(), res.dimensions(), mps1_ref.data(), mps1_ref.dimensions(), mps2_ref.data(), mps2_ref.dimensions(), idx);
        return res;
    }

    template<typename mps_type>
    [[nodiscard]] Eigen::Tensor<typename mps_type::Scalar, 2> contract_mps_partial(const TensorRead<mps_type> &mps, std::array<long, 2> idx) {
        static_assert(mps_type::NumIndices == 3 and "Wrong mps tensor rank != 3");
        const auto         &mps_ref = static_cast<const mps_type &>(mps);
        std::array<long, 2> dims    = {};
        if(idx[1] == 1)
            dims = {mps_ref.dimension(2), mps_ref.dimension(2)};
        else if(idx[1] == 2)
            dims = {mps_ref.dimension(1), mps_ref.dimension(1)};
        else
            throw except::runtime_error("Expected idx [0,1] or [0,2]. Got: {}", idx);
        Eigen::Tensor<typename mps_type::Scalar, 2> res(dims);
        contract_mps_mps_partial(res.data(), res.dimensions(), mps_ref.data(), mps_ref.dimensions(), mps_ref.data(), mps_ref.dimensions(), idx);
        return res;
    }

}