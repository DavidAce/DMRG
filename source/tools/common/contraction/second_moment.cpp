// template<typename mps_type, typename mpo_type, typename env_type>
// extern double tools::common::contraction::second_moment(const mps_type &mps, const mpo_type &mpo, const env_type &envL, const env_type &envR) {
//    // This measures the second moment <M²> of some multisite operator M given multisite mps', mpos and corresponding environments.
//    // This is usually the second moment of the hamiltonian <H²>
//    // Note that the environments must contain the correct type of mpos
//    using Scalar    = typename mps_type::Scalar;
//    double log2chiL = std::log2(mps.dimension(1));
//    double log2chiR = std::log2(mps.dimension(2));
//    double log2spin = std::log2(mps.dimension(0));
//    if(mps.dimension(1) != envL.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions()));
//    if(mps.dimension(2) != envR.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions()));
//    if(mps.dimension(0) != mpo.dimension(2)) throw std::runtime_error(fmt::format("Dimension mismatch mps {} and mpo {}", mps.dimensions(),
//    mpo.dimensions()));
//
//    Eigen::Tensor<Scalar, 0> M2;
//    /* clang-format off */
//    if(log2spin >= std::max(log2chiL, log2chiR)) {
//        if(log2chiL > log2chiR) {
//            //            tools::log->trace("H2 path: log2spin > std::max(log2chiL , log2chiR)  and  log2chiL > log2chiR ");
//            Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(tenx::array3{1, 0, 2});
//            M2.device(tenx::omp::getDevice()) =
//                mps_shuffled
//                    .contract(envL, tenx::idx({0}, {0}))
//                    .contract(mpo,  tenx::idx({0, 3}, {2, 0}))
//                    .contract(envR, tenx::idx({0, 3}, {0, 2}))
//                    .contract(mpo , tenx::idx({2, 1, 4}, {2, 0, 1}))
//                    .contract(mps_shuffled.conjugate(), tenx::idx({2, 0, 1}, {1, 0, 2}));
//        }
//        else {
//            //            tools::log->trace("H2 path: log2spin >= std::max(log2chiL , log2chiR) and  log2chiL <= log2chiR ");
//            Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(tenx::array3{2, 0, 1});
//            M2.device(tenx::omp::getDevice()) =
//                mps_shuffled
//                    .contract(envR, tenx::idx({0}, {0}))
//                    .contract(mpo,  tenx::idx({0, 3}, {2, 1}))
//                    .contract(envL, tenx::idx({0, 3}, {0, 2}))
//                    .contract(mpo,  tenx::idx({2, 4, 1}, {2, 0, 1}))
//                    .contract(mps_shuffled.conjugate(), tenx::idx({2, 1, 0}, {1, 2, 0}));
//        }
//    } else {
//        //        tools::log->trace("H2 path: log2spin < std::max(log2chiL , log2chiR)");
//        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(tenx::array3{1, 0, 2});
//        M2.device(tenx::omp::getDevice()) =
//            mps_shuffled.contract(envL, tenx::idx({0}, {0}))
//                .contract(mpo,  tenx::idx({0, 3}, {2, 0}))
//                .contract(mpo,  tenx::idx({4, 2}, {2, 0}))
//                .contract(envR, tenx::idx({0, 2, 3}, {0, 2, 3}))
//                .contract(mps_shuffled.conjugate(), tenx::idx({1, 0, 2}, {1, 0, 2}));
//    }
//    /* clang-format on */
//    if(abs(std::imag(M2(0))) > 1e-10) {
//        throw std::runtime_error(fmt::format("Second moment has an imaginary part: {:.16f} + i {:.16f}", std::real(M2(0)), std::imag(M2(0))));
//    }
//    double moment = std::real(M2(0));
//    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("Second moment is invalid: {}", moment));
//    return moment;
//}

// template double tools::common::moments::second_moment(const T3<real> &, const T4<real> &, const T4<real> &, const T4<real> &);
// template double tools::common::moments::second_moment(const T3<cplx> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);
// template double tools::common::moments::second_moment(const TM<T3<real>> &, const T4<real> &, const T4<real> &, const T4<real> &);
// template double tools::common::moments::second_moment(const TM<T3<cplx>> &, const T4<cplx> &, const T4<cplx> &, const T4<cplx> &);