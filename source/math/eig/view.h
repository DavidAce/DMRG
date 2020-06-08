#pragma once

#include "solution.h"
#include <general/nmspc_tensor_extra.h>

namespace eig::view {

    template<typename T>
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename T>
    using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;


    template<typename Scalar>
    Eigen::Map<VectorType<Scalar>> get_eigvals(eig::solution &result) {
        if(not result.meta.eigvals_found)
            throw std::runtime_error(fmt::format("Results have not been obtained yet: eigvals found [{}]", result.meta.eigvals_found));
        if constexpr(std::is_same_v<Scalar, real>)
            if(not result.eigvals_are_real()) throw std::runtime_error("Can't view real eigenvalues: solution has complex eigenvalues");

        auto & eigvals = result.get_eigvals<Scalar>();
        if(eigvals.empty()) throw std::runtime_error("The requested eigenvalues are empty. Did you request the correct type?");
        return Eigen::Map<VectorType<Scalar>>(eigvals.data(), static_cast<long>(eigvals.size()));
    }

    template<typename Scalar, typename integral_type = long, typename = std::enable_if<std::is_integral_v<integral_type>>>
    Scalar get_eigval(eig::solution &result, integral_type num = 0) {
        return get_eigvals<Scalar>(result)(static_cast<long>(num));
    }


    template<typename Scalar>
    Eigen::Map<MatrixType<Scalar>> get_eigvecs(eig::solution &result, Side side = Side::R) {
        if(side == Side::R and not result.meta.eigvecsR_found)
            throw std::runtime_error(fmt::format("Results have not been obtained yet: eigvecs R found [{}]", result.meta.eigvecsR_found));
        else if(side == Side::L and not result.meta.eigvecsL_found)
            throw std::runtime_error(fmt::format("Results have not been obtained yet: eigvecs L found [{}]", result.meta.eigvecsL_found));

        if(side == Side::LR) throw std::runtime_error("Cannot get both L and R eigenvectors. Choose one");

        if constexpr(std::is_same_v<Scalar, real>)
            if(not result.eigvecs_are_real()) throw std::runtime_error("Can't view real eigenvectors: solution has complex eigenvectors");
        auto & eigvecs = result.get_eigvecs<Scalar>(side);
        if(eigvecs.empty()) throw std::runtime_error("The requested eigenvectors are empty. Did you request the correct type?");
        auto rows = result.meta.rows;
        auto cols = result.meta.cols;
        if(rows * cols != static_cast<eig::size_type>(eigvecs.size())) throw std::logic_error("Size mismatch in results");
        return Eigen::Map<MatrixType<Scalar>>(eigvecs.data(), result.meta.rows, result.meta.cols);
    }

    template<typename Scalar, typename integral_type = long, typename = std::enable_if<std::is_integral_v<integral_type>>>
    Eigen::Map<VectorType<Scalar>> get_eigvec(eig::solution &result, integral_type num = 0, Side side = Side::R) {
        auto eigvecmap = get_eigvecs<Scalar>(result,side).col(static_cast<long>(num));
        return Eigen::Map<VectorType<Scalar>>(eigvecmap.data(), eigvecmap.size());
    }

    template<typename Scalar, typename integral_type = long, typename = std::enable_if<std::is_integral_v<integral_type>>>
    Eigen::TensorMap<Eigen::Tensor<Scalar, 3>> get_eigvec(eig::solution &result, const Eigen::DSizes<long, 3> &dims, integral_type num = 0, Side side = Side::R) {
        if(result.meta.rows != dims[0] * dims[1] * dims[2])
            throw std::range_error(fmt::format("Given tensor dimensions do not match eigenvector size: size {} != {} * {} * {} = {}", result.meta.rows, dims[0],
                                               dims[1], dims[3], dims[0] * dims[1] * dims[2]));
        auto eigvecmap = get_eigvec<Scalar>(result, num, side);
        return Eigen::TensorMap<Eigen::Tensor<Scalar, 3>>(eigvecmap.data(), dims);
    }




}
