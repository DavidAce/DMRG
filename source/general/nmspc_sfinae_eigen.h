#pragma once
#include <Eigen/Core>
#include <general/eigen_tensor_fwd_decl.h>
#include <general/nmspc_sfinae.h>

/*!
 * \brief A collection of type-detection and type-analysis utilities using SFINAE on Eigen types
 */
namespace sfinae::eigen {

    template<typename T>
    using is_eigen_matrix = std::is_base_of<Eigen::MatrixBase<std::decay_t<T>>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_matrix_v = is_eigen_matrix<T>::value;
    template<typename T>
    using is_eigen_array = std::is_base_of<Eigen::ArrayBase<std::decay_t<T>>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_array_v = is_eigen_array<T>::value;
    template<typename T>
    using is_eigen_tensor = std::is_base_of<Eigen::TensorBase<std::decay_t<T>, Eigen::ReadOnlyAccessors>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_tensor_v = is_eigen_tensor<T>::value;
    template<typename T>
    using is_eigen_dense = std::is_base_of<Eigen::DenseBase<std::decay_t<T>>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_dense_v = is_eigen_dense<T>::value;
    template<typename T>
    using is_eigen_map = std::is_base_of<Eigen::MapBase<std::decay_t<T>, Eigen::ReadOnlyAccessors>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_map_v = is_eigen_map<T>::value;
    template<typename T>
    using is_eigen_plain = std::is_base_of<Eigen::PlainObjectBase<std::decay_t<T>>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_plain_v = is_eigen_plain<T>::value;
    template<typename T>
    using is_eigen_base = std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>;
    template<typename T>
    inline constexpr bool is_eigen_base_v = is_eigen_base<T>::value;
    template<typename T>
    struct is_eigen_core : public std::false_type {};
    template<typename T, int rows, int cols, int StorageOrder>
    struct is_eigen_core<Eigen::Matrix<T, rows, cols, StorageOrder>> : public std::true_type {};
    template<typename T, int rows, int cols, int StorageOrder>
    struct is_eigen_core<Eigen::Array<T, rows, cols, StorageOrder>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_eigen_core_v = is_eigen_core<T>::value;
    template<typename T>
    struct is_eigen_any {
        static constexpr bool value = is_eigen_base<T>::value or is_eigen_tensor<T>::value;
    };
    template<typename T>
    inline constexpr bool is_eigen_any_v = is_eigen_any<T>::value;
    template<typename T>
    struct is_eigen_contiguous {
        static constexpr bool value = is_eigen_any<T>::value and has_data<T>::value;
    };
    template<typename T>
    inline constexpr bool is_eigen_contiguous_v = is_eigen_contiguous<T>::value;
    template<typename T>
    class is_eigen_1d {
        private:
        template<typename U>
        static constexpr auto test() {
            if constexpr(is_eigen_map<U>::value) return test<typename U::PlainObject>();
            if constexpr(is_eigen_dense<U>::value) return U::RowsAtCompileTime == 1 or U::ColsAtCompileTime == 1;
            if constexpr(is_eigen_tensor<U>::value and has_NumIndices<U>::value)
                return U::NumIndices == 1;
            else
                return false;
        }

        public:
        static constexpr bool value = test<T>();
    };
    template<typename T>
    inline constexpr bool is_eigen_1d_v = is_eigen_1d<T>::value;

    template<typename T>
    class is_eigen_colmajor {
        template<typename U>
        static constexpr bool test() {
            if constexpr(is_eigen_base<U>::value) return not U::IsRowMajor;
            if constexpr(is_eigen_tensor<U>::value)
                return Eigen::ColMajor == static_cast<Eigen::StorageOptions>(U::Layout);
            else
                return false;
        }

        public:
        static constexpr bool value = test<T>();
    };
    template<typename T>
    inline constexpr bool is_eigen_colmajor_v = is_eigen_colmajor<T>::value;

    template<typename T>
    class is_eigen_rowmajor {
        template<typename U>
        static constexpr bool test() {
            if constexpr(is_eigen_base<U>::value) return U::IsRowMajor;
            if constexpr(is_eigen_tensor<U>::value)
                return Eigen::RowMajor == static_cast<Eigen::StorageOptions>(U::Layout);
            else
                return false;
        }

        public:
        static constexpr bool value = test<T>();
    };
    template<typename T>
    inline constexpr bool is_eigen_rowmajor_v = is_eigen_rowmajor<T>::value;

}
