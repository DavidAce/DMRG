#pragma once

namespace tools::common::moments{
    template<typename mps_type, typename mpo_type, typename env_type>
    double first(const mps_type &mps,const mpo_type &mpo, const env_type &envL, const env_type &envR);
    template<typename mps_type, typename mpo_type, typename env_type>
    double second(const mps_type &mps,const mpo_type &mpo, const env_type &envL, const env_type &envR);

    // Extern templates
    using cplx = std::complex<double>;
    using real = double;
    template<typename Scalar>
    using T3 = Eigen::Tensor<Scalar, 3>;
    template<typename Scalar>
    using T4 = Eigen::Tensor<Scalar, 4>;
    template<typename T>
    using TM = Eigen::TensorMap<T>;
//    extern template double tools::common::moments::first(const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
    extern template double tools::common::moments::first(const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
//    extern template double tools::common::moments::first(const TM<T3<real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
//    extern template double tools::common::moments::first(const TM<T3<cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
//    extern template double tools::common::moments::second(const T3<real> &, const T4<real> &, const T4<real> &, const T4<real> &);
    extern template double tools::common::moments::second(const T3<cplx> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);
//    extern template double tools::common::moments::second(const TM<T3<real>> &, const T4<real> &, const T4<real> &, const T4<real> &);
//    extern template double tools::common::moments::second(const TM<T3<cplx>> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);


}