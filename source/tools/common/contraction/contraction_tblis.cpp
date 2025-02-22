#include "../contraction.h"
#include "io/fmt_custom.h"
#include "math/tenx.h"
#include "tid/tid.h"
#include <fmt/ranges.h>
#include <omp.h>
#if defined(DMRG_ENABLE_TBLIS)
    #include <tblis/tblis.h>
    #include <tblis/util/thread.h>
    #include <tci/tci_config.h>
#endif
#if defined(DMRG_SAVE_CONTRACTION)
    #include <h5pp/h5pp.h>
#endif
#if defined(DMRG_BENCH_CONTRACTION)
    #include "math/num.h"
    #include "tid/tid.h"
namespace settings {
    [[maybe_unused]] constexpr static bool bench_expval = true;
}
#else
namespace settings {
    [[maybe_unused]] constexpr static bool bench_expval = false;
}
#endif

using namespace tools::common::contraction;

/* clang-format off */
template<typename Scalar>
Scalar tools::common::contraction::expectation_value(const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                                     const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                                                     const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                                                     const Scalar * const envR_ptr, std::array<long,3> envR_dims){

    std::string bench_suffix;
#if defined(DMRG_BENCH_CONTRACTION)
    if constexpr (settings::bench_expval){
        auto spin = mps_dims[0];
        auto chiL = num::next_multiple<long>(mps_dims[1], 5l);
        auto chiR = num::next_multiple<long>(mps_dims[2], 5l);
        auto mdim = mpo_dims[0];
        bench_suffix = fmt::format("_mps-[{},{},{}]_mpo-[{}]", spin, chiL, chiR, mdim);
    }
#endif
    auto t_expval = tid::tic_token(fmt::format("expval{}",bench_suffix), tid::level::highest);

    // This measures the expectation value of some multisite mps with respect to some mpo operator and corresponding environments.
    // This is usually the energy E = <psi|H|psi> or variance V = <psi|(H-E)²|psi>
    // Note that the environments must contain the correct type of mpos
    auto mps = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(mps_ptr,mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<Scalar,4>>(mpo_ptr,mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envR_ptr,envR_dims);

    assert(mps.dimension(1)  == envL.dimension(0));
    assert(mps.dimension(2)  == envR.dimension(0));
    assert(mps.dimension(0)  == mpo.dimension(2));
    assert(envL.dimension(2) == mpo.dimension(0));
    assert(envR.dimension(2) == mpo.dimension(1));

    Eigen::Tensor<Scalar, 0> expval;
    auto &threads = tenx::threads::get();
    expval.device(*threads->dev) =
        envL
            .contract(mps,             tenx::idx({0}, {1}))
            .contract(mpo,             tenx::idx({2, 1}, {2, 0}))
            .contract(mps.conjugate(), tenx::idx({3, 0}, {0, 1}))
            .contract(envR,            tenx::idx({0, 2, 1}, {0, 1, 2}));

    Scalar result = expval(0);
    #if defined(DMRG_SAVE_CONTRACTION)
    {
        t_expval.toc();
        auto file       = h5pp::File("dmrg-contractions.h5", h5pp::FilePermission::READWRITE);
        auto group_num  = 0;
        auto group_name = fmt::format("contraction_{}", group_num);
        while(file.linkExists(group_name)) group_name = fmt::format("contraction_{}", ++group_num);
        file.writeDataset(mps, fmt::format("{}/mps", group_name));
        file.writeDataset(mpo, fmt::format("{}/mpo", group_name));
        file.writeDataset(envL, fmt::format("{}/envL", group_name));
        file.writeDataset(envR, fmt::format("{}/envR", group_name));
    }
    #endif

    return result;
}

template fp64 tools::common::contraction::expectation_value(const fp64 * const mps_ptr,  std::array<long,3> mps_dims,
                                                            const fp64 * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                            const fp64 * const envL_ptr, std::array<long,3> envL_dims,
                                                            const fp64 * const envR_ptr, std::array<long,3> envR_dims);
template cx64 tools::common::contraction::expectation_value(const cx64 * const mps_ptr,  std::array<long,3> mps_dims,
                                                            const cx64 * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                            const cx64 * const envL_ptr, std::array<long,3> envL_dims,
                                                            const cx64 * const envR_ptr, std::array<long,3> envR_dims);
/* clang-format on */

/* clang-format off */
template<typename Scalar>
Scalar tools::common::contraction::expectation_value(const Scalar * const bra_ptr, std::array<long,3> bra_dims,
                                                     const Scalar * const ket_ptr, std::array<long,3> ket_dims,
                                                     const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                                                     const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                                                     const Scalar * const envR_ptr, std::array<long,3> envR_dims){

    std::string bench_suffix;
#if defined(DMRG_BENCH_CONTRACTION)
    if constexpr (settings::bench_expval){
        auto spin = mps_dims[0];
        auto chiL = num::next_multiple<long>(mps_dims[1], 5l);
        auto chiR = num::next_multiple<long>(mps_dims[2], 5l);
        auto mdim = mpo_dims[0];
        bench_suffix = fmt::format("_mps-[{},{},{}]_mpo-[{}]", spin, chiL, chiR, mdim);
    }
#endif
    auto t_expval = tid::tic_token(fmt::format("expval{}",bench_suffix), tid::level::highest);

    // This measures the expectation value of some multisite mps with respect to some mpo operator and corresponding environments.
    // This is usually the energy E = <psi|H|psi> or variance V = <psi|(H-E)²|psi>
    // Note that the environments must contain the correct type of mpos
    auto bra = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(bra_ptr,bra_dims);
    auto ket = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(ket_ptr,ket_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<Scalar,4>>(mpo_ptr,mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envR_ptr,envR_dims);

    assert(bra.dimension(1)  == envL.dimension(1));
    assert(ket.dimension(1)  == envL.dimension(0));
    assert(bra.dimension(2)  == envR.dimension(1));
    assert(ket.dimension(2)  == envR.dimension(0));
    assert(bra.dimension(0)  == mpo.dimension(3));
    assert(ket.dimension(0)  == mpo.dimension(2));
    assert(envL.dimension(2) == mpo.dimension(0));
    assert(envR.dimension(2) == mpo.dimension(1));

    Eigen::Tensor<Scalar, 0> expval;
    auto &threads = tenx::threads::get();
    expval.device(*threads->dev) =
        envL
            .contract(ket,             tenx::idx({0}, {1}))
            .contract(mpo,             tenx::idx({2, 1}, {2, 0}))
            .contract(bra.conjugate(), tenx::idx({3, 0}, {0, 1}))
            .contract(envR,            tenx::idx({0, 2, 1}, {0, 1, 2}));

    Scalar result = expval(0);
    #if defined(DMRG_SAVE_CONTRACTION)
    {
        t_expval.toc();
        auto file       = h5pp::File("dmrg-contractions.h5", h5pp::FilePermission::READWRITE);
        auto group_num  = 0;
        auto group_name = fmt::format("contraction_{}", group_num);
        while(file.linkExists(group_name)) group_name = fmt::format("contraction_{}", ++group_num);
        file.writeDataset(mps, fmt::format("{}/mps", group_name));
        file.writeDataset(mpo, fmt::format("{}/mpo", group_name));
        file.writeDataset(envL, fmt::format("{}/envL", group_name));
        file.writeDataset(envR, fmt::format("{}/envR", group_name));
    }
    #endif

    return result;
}

template fp64 tools::common::contraction::expectation_value(const fp64 * const bra_ptr,  std::array<long,3> bra_dims,
                                                            const fp64 * const ket_ptr,  std::array<long,3> ket_dims,
                                                            const fp64 * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                            const fp64 * const envL_ptr, std::array<long,3> envL_dims,
                                                            const fp64 * const envR_ptr, std::array<long,3> envR_dims);
template cx64 tools::common::contraction::expectation_value(const cx64 * const bra_ptr,  std::array<long,3> bra_dims,
                                                            const cx64 * const ket_ptr,  std::array<long,3> ket_dims,
                                                            const cx64 * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                            const cx64 * const envL_ptr, std::array<long,3> envL_dims,
                                                            const cx64 * const envR_ptr, std::array<long,3> envR_dims);


#if defined(DMRG_ENABLE_TBLIS)
std::string get_tblis_arch() {
    #if defined(__GNUC__)
    if(__builtin_cpu_supports("x86-64-v4")) return "skx";
    if(__builtin_cpu_is("znver3") or __builtin_cpu_is("znver2") or __builtin_cpu_is("znver1")) return "zen";
    if(__builtin_cpu_supports("x86-64-v3")) return "haswell";
    #endif
    return "haswell";
}

template<typename ea_type, typename eb_type, typename ec_type>
void contract_tblis(const TensorRead<ea_type> &ea, const TensorRead<eb_type> &eb, TensorWrite<ec_type> &ec, const tblis::label_vector &la,
                    const tblis::label_vector &lb, const tblis::label_vector &lc, const tblis::tblis_config_s *tblis_config) {
    const auto &ea_ref = static_cast<const ea_type &>(ea);
    const auto &eb_ref = static_cast<const eb_type &>(eb);
    auto       &ec_ref = static_cast<ec_type &>(ec);

    tblis::len_vector da, db, dc;
    da.assign(ea_ref.dimensions().begin(), ea_ref.dimensions().end());
    db.assign(eb_ref.dimensions().begin(), eb_ref.dimensions().end());
    dc.assign(ec_ref.dimensions().begin(), ec_ref.dimensions().end());

    auto                     ta    = tblis::varray_view<const typename ea_type::Scalar>(da, ea_ref.data(), tblis::COLUMN_MAJOR);
    auto                     tb    = tblis::varray_view<const typename eb_type::Scalar>(db, eb_ref.data(), tblis::COLUMN_MAJOR);
    auto                     tc    = tblis::varray_view<typename ec_type::Scalar>(dc, ec_ref.data(), tblis::COLUMN_MAJOR);
    typename ea_type::Scalar alpha = 1.0;
    typename ec_type::Scalar beta  = 0.0;

    tblis::tblis_tensor A_s(alpha, ta);
    tblis::tblis_tensor B_s(tb);
    tblis::tblis_tensor C_s(beta, tc);

    tblis_tensor_mult(nullptr, tblis_config, &A_s, la.c_str(), &B_s, lb.c_str(), &C_s, lc.c_str());
}
#endif

/* clang-format off */
template<typename Scalar>
void tools::common::contraction::matrix_vector_product(      Scalar * res_ptr,
                                                       const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                                       const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                                                       const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                                                       const Scalar * const envR_ptr, std::array<long,3> envR_dims){

//    auto t_matvec = tid::tic_token("matrix_vector_product", tid::level::extra);

    // This applies the mpo's with corresponding environments to local multisite mps
    // This is usually the operation H|psi>  or H²|psi>
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,mps_dims);
    auto mps = Eigen::TensorMap<const  Eigen::Tensor<Scalar,3>>(mps_ptr,mps_dims);
    auto mpo = Eigen::TensorMap<const  Eigen::Tensor<Scalar,4>>(mpo_ptr,mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envR_ptr,envR_dims);

    if(mps.dimension(1) != envL.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions());
    if(mps.dimension(2) != envR.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions());
    if(mps.dimension(0) != mpo.dimension(2))  throw except::runtime_error("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions());
    if(envL.dimension(2) != mpo.dimension(0)) throw except::runtime_error("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions());
    if(envR.dimension(2) != mpo.dimension(1)) throw except::runtime_error("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions());

#if defined(DMRG_ENABLE_TBLIS)
    if constexpr(std::is_same_v<Scalar, fp64>){
        static const tblis::tblis_config_s *tblis_config = tblis::tblis_get_config(get_tblis_arch().data());
        #if defined(TCI_USE_OPENMP_THREADS) && defined(_OPENMP)
        tblis_set_num_threads(static_cast<unsigned int>(omp_get_max_threads()));
        #endif
        if (mps.dimension(1) >= mps.dimension(2)){
            Eigen::Tensor<Scalar, 4> mpsenvL(mps.dimension(0), mps.dimension(2), envL.dimension(1), envL.dimension(2));
            // Eigen::Tensor<Scalar, 4> mpsenvLmpo(mps.dimension(2), envL.dimension(1), mpo.dimension(1), mpo.dimension(3));
            Eigen::Tensor<Scalar, 4> mpsenvLmpo_alt(mpo.dimension(1), mpo.dimension(3), mps.dimension(2), envL.dimension(1));

            contract_tblis(mps, envL, mpsenvL, "afb", "fcd", "abcd", tblis_config);
            // contract_tblis(mpsenvL, mpo, mpsenvLmpo, "qijr", "rkql", "ijkl", tblis_config);
            // contract_tblis(mpsenvLmpo, envR, res, "qjri", "qkr", "ijk", tblis_config);
            contract_tblis(mpo, mpsenvL, mpsenvLmpo_alt, "qhri", "rgjq", "higj", tblis_config);
            contract_tblis(mpsenvLmpo_alt, envR, res, "higj", "gkh", "ijk", tblis_config);
        }
        else{
            Eigen::Tensor<Scalar, 4> mpsenvR(mps.dimension(0), mps.dimension(1), envR.dimension(1), envR.dimension(2));
            Eigen::Tensor<Scalar, 4> mpsenvRmpo(mps.dimension(1), envR.dimension(1), mpo.dimension(0), mpo.dimension(3));
            contract_tblis(mps, envR, mpsenvR, "abf", "fcd", "abcd", tblis_config);
            contract_tblis(mpsenvR, mpo, mpsenvRmpo, "qijk", "rkql", "ijrl", tblis_config);
            contract_tblis(mpsenvRmpo, envL, res, "qkri", "qjr", "ijk", tblis_config);
        }
    }else{
        auto &threads = tenx::threads::get();
        if (mps.dimension(1) >= mps.dimension(2)){
            res.device(*threads->dev) = mps
                                     .contract(envL, tenx::idx({1}, {0}))
                                     .contract(mpo,  tenx::idx({3, 0}, {0, 2}))
                                     .contract(envR, tenx::idx({0, 2}, {0, 2}))
                                     .shuffle(tenx::array3{1, 0, 2});
        }else{
            res.device(*threads->dev) = mps
                                     .contract(envR, tenx::idx({2}, {0}))
                                     .contract(mpo,  tenx::idx({3, 0}, {1, 2}))
                                     .contract(envL, tenx::idx({0, 2}, {0, 2}))
                                     .shuffle(tenx::array3{1, 2, 0});
        }
    }
    #else
    auto &threads = tenx::threads::get();
    if (mps.dimension(1) >= mps.dimension(2)){
        res.device(*threads->dev) = mps
                                 .contract(envL, tenx::idx({1}, {0}))
                                 .contract(mpo,  tenx::idx({3, 0}, {0, 2}))
                                 .contract(envR, tenx::idx({0, 2}, {0, 2}))
                                 .shuffle(tenx::array3{1, 0, 2});
    }else{
        res.device(*threads->dev) = mps
                                 .contract(envR, tenx::idx({2}, {0}))
                                 .contract(mpo,  tenx::idx({3, 0}, {1, 2}))
                                 .contract(envL, tenx::idx({0, 2}, {0, 2}))
                                 .shuffle(tenx::array3{1, 2, 0});
    }
    #endif
}

using namespace tools::common::contraction;
template void tools::common::contraction::matrix_vector_product(      cx64 *       res_ptr,
                                                                const cx64 * const mps_ptr, std::array<long,3> mps_dims,
                                                                const cx64 * const mpo_ptr, std::array<long,4> mpo_dims,
                                                                const cx64 * const envL_ptr, std::array<long,3> envL_dims,
                                                                const cx64 * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_vector_product(      fp64 *       res_ptr,
                                                                const fp64 * const mps_ptr, std::array<long,3> mps_dims,
                                                                const fp64 * const mpo_ptr, std::array<long,4> mpo_dims,
                                                                const fp64 * const envL_ptr, std::array<long,3> envL_dims,
                                                                const fp64 * const envR_ptr, std::array<long,3> envR_dims);


template<typename Scalar, typename mpo_type>
void tools::common::contraction::matrix_vector_product(Scalar * res_ptr,
                           const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                           const std::vector<mpo_type> & mpos_shf,
                           const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                           const Scalar * const envR_ptr, std::array<long,3> envR_dims) {
    // Make sure the mpos are pre-shuffled. If not, shuffle and call this function again
    bool is_shuffled = mpos_shf.front().dimension(2) == envL_dims[2] and mpos_shf.back().dimension(3) == envR_dims[2];
    if(not is_shuffled){
        // mpos_shf are not actually shuffled. Let's shuffle.
        std::vector<Eigen::Tensor<Scalar, 4>> mpos_really_shuffled;
        for (const auto & mpo : mpos_shf) {
            mpos_really_shuffled.emplace_back(mpo.shuffle(tenx::array4{2, 3, 0, 1}));
        }
        return matrix_vector_product(res_ptr, mps_ptr, mps_dims, mpos_really_shuffled, envL_ptr, envL_dims, envR_ptr, envR_dims);
    }


    auto &threads  = tenx::threads::get();
    auto mps_out = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,mps_dims);
    auto mps_in  = Eigen::TensorMap<const  Eigen::Tensor<Scalar,3>>(mps_ptr,mps_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(envR_ptr,envR_dims);

    if(mps_in.dimension(1) != envL.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envL {}", mps_in.dimensions(), envL.dimensions());
    if(mps_in.dimension(2) != envR.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envR {}", mps_in.dimensions(), envR.dimensions());


    auto  L        = mpos_shf.size();

    auto mpodimprod = [&](size_t fr, size_t to) -> long {
        long prod = 1;
        if(fr == -1ul) fr = 0;
        if(to == 0 or to == -1ul) return prod;
        for(size_t idx = fr; idx < to; ++idx) {
            if(idx >= mpos_shf.size()) break;
            prod *= mpos_shf[idx].dimension(1);
        }
        return prod;
    };

    // At best, the number of operations for contracting left-to-right and right-to-left are equal.
    // Since the site indices are contracted left to right, we do not need any shuffles in this direction.

    // Contract left to right
    #if defined(DMRG_ENABLE_TBLIS)
    static const tblis::tblis_config_s *tblis_config = tblis::tblis_get_config(get_tblis_arch().data());
    #if defined(TCI_USE_OPENMP_THREADS) && defined(_OPENMP)
    tblis_set_num_threads(static_cast<unsigned int>(omp_get_max_threads()));
    // if(omp_get_max_active_levels() <= 1) omp_set_max_active_levels(2);
    // printf("tblis: %d | omp_get_num_threads: %d / %d | omp_get_active_level %d | omp_in_parallel %d\n",
            // tblis_get_num_threads(), omp_get_num_threads(), omp_get_max_threads(), omp_get_active_level(), omp_in_parallel());
    #endif
    #endif
    auto d0       = mpodimprod(0, 1); // Split 0 --> 0,1
    auto d1       = mpodimprod(1, L); // Split 0 --> 0,1
    auto d2       = mps_in.dimension(2);
    auto d3       = envL.dimension(1);
    auto d4       = envL.dimension(2);
    auto d5       = 1l; // A new dummy index
    auto new_shp6 = tenx::array6{d0, d1, d2, d3, d4, d5};
    auto  mps_tmp1 = Eigen::Tensor<Scalar, 6>();
    auto  mps_tmp2 = Eigen::Tensor<Scalar, 6>();
    mps_tmp1.resize(tenx::array6{d0, d1, d2, d3, d5, d4});
    #if defined(DMRG_ENABLE_TBLIS)
    if constexpr(std::is_same_v<Scalar, fp64>) {
        auto mps_tmp1_map4 = Eigen::TensorMap<Eigen::Tensor<Scalar, 4>>(mps_tmp1.data(), std::array{d0 * d1, d2, d3, d4 * d5});
        contract_tblis(mps_in, envL, mps_tmp1_map4, "afb", "fcd", "abcd", tblis_config);
    } else
    #endif
    {
        mps_tmp1.device(*threads->dev) = mps_in.contract(envL, tenx::idx({1}, {0})).reshape(new_shp6).shuffle(tenx::array6{0, 1, 2, 3, 5, 4});
    }
    for(size_t idx = 0; idx < L; ++idx) {
        const auto &mpo = mpos_shf[idx];
        // Set up the dimensions for the reshape after the contraction
        d0       = mpodimprod(idx + 1, idx + 2); // if idx == k, this has the mpo at idx == k+1
        d1       = mpodimprod(idx + 2, L);       // if idx == 0,  this has the mpos at idx == k+2...L-1
        d2       = mps_tmp1.dimension(2);
        d3       = mps_tmp1.dimension(3);
        d4       = mpodimprod(0, idx + 1); // if idx == 0, this has the mpos at idx == 0...k (i.e. including the one from the current iteration)
        d5       = mpo.dimension(3);       // The virtual bond of the current mpo
        #if defined(DMRG_ENABLE_TBLIS)
        if constexpr(std::is_same_v<Scalar, fp64>) {
            auto md  = mps_tmp1.dimensions();
            new_shp6 = tenx::array6{d0, d1, d2, d3, d4, d5};
            mps_tmp2.resize(new_shp6);
            auto map_shp6 = tenx::array6{md[1], md[2], md[3], md[4], mpo.dimension(1), mpo.dimension(3)};
            auto mps_tmp2_map = Eigen::TensorMap<Eigen::Tensor<Scalar,6>>(mps_tmp2.data(), map_shp6);
            contract_tblis(mps_tmp1, mpo, mps_tmp2_map, "qbcder", "qfrg", "bcdefg", tblis_config);
            mps_tmp1 = std::move(mps_tmp2);
        } else
        #endif
        {
            new_shp6 = tenx::array6{d0, d1, d2, d3, d4, d5};
            mps_tmp2.resize(new_shp6);
            mps_tmp2.device(*threads->dev) = mps_tmp1.contract(mpo, tenx::idx({0, 5}, {0, 2})).reshape(new_shp6);
            mps_tmp1                       = std::move(mps_tmp2);
        }
    }
    d0 = mps_tmp1.dimension(0) * mps_tmp1.dimension(1) * mps_tmp1.dimension(2); // idx 0 and 1 should have dim == 1
    d1 = mps_tmp1.dimension(3);
    d2 = mps_tmp1.dimension(4);
    d3 = mps_tmp1.dimension(5);
    #if defined(DMRG_ENABLE_TBLIS)
    if constexpr(std::is_same_v<Scalar, fp64>) {
        auto mps_tmp1_map4 = Eigen::TensorMap<Eigen::Tensor<Scalar, 4>>(mps_tmp1.data(), std::array{d0, d1, d2, d3});
        contract_tblis(mps_tmp1_map4, envR, mps_out, "qjir", "qkr", "ijk", tblis_config);
    } else
    #endif
    {
        mps_out.device(*threads->dev) = mps_tmp1.reshape(tenx::array4{d0, d1, d2, d3}).contract(envR, tenx::idx({0, 3}, {0, 2})).shuffle(tenx::array3{1, 0, 2});
    }
}

using namespace tools::common::contraction;
template void tools::common::contraction::matrix_vector_product(      cx64 *       res_ptr,
                                                                const cx64 * const mps_ptr, std::array<long,3> mps_dims,
                                                                const std::vector<Eigen::Tensor<cx64, 4>> & mpos_shf,
                                                                const cx64 * const envL_ptr, std::array<long,3> envL_dims,
                                                                const cx64 * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_vector_product(      fp64 *       res_ptr,
                                                                const fp64 * const mps_ptr, std::array<long,3> mps_dims,
                                                                const std::vector<Eigen::Tensor<fp64, 4>> & mpos_shf,
                                                                const fp64 * const envL_ptr, std::array<long,3> envL_dims,
                                                                const fp64 * const envR_ptr, std::array<long,3> envR_dims);

template<typename Scalar>
void  tools::common::contraction::contract_mps_bnd(      Scalar * res_ptr      , std::array<long,3> res_dims,
                                                   const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                                   const Scalar * const bnd_ptr, std::array<long,1> bnd_dims){
//    auto t_con = tid::tic_token("contract_mps_bnd", tid::level::highest);
    if(mps_dims[2] != bnd_dims[0]) throw except::runtime_error("Dimension mismatch mps {} (idx 2) and bnd {} (idx 0)", mps_dims, bnd_dims);
    if(mps_dims != res_dims) throw except::runtime_error("Dimension mismatch mps {} and res {}", mps_dims, res_dims);

    auto res_mat = Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic>>(res_ptr, res_dims[0] * res_dims[1], res_dims[2]);
    auto mps_mat = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic>>(mps_ptr, mps_dims[0] * mps_dims[1], mps_dims[2]);
    auto bnd_mat = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic, 1>>(bnd_ptr, bnd_dims[0]);
    res_mat.noalias() = mps_mat * bnd_mat.asDiagonal(); // calls gemm
}
template void tools::common::contraction::contract_mps_bnd(      cx64 *       res_ptr, std::array<long,3> res_dims,
                                                           const cx64 * const mps_ptr, std::array<long,3> mps_dims,
                                                           const cx64 * const bnd_ptr, std::array<long,1> bnd_dims);

template void tools::common::contraction::contract_mps_bnd(      fp64 *       res_ptr, std::array<long,3> res_dims,
                                                           const fp64 * const mps_ptr, std::array<long,3> mps_dims,
                                                           const fp64 * const bnd_ptr, std::array<long,1> bnd_dims);


template<typename Scalar>
void  tools::common::contraction::contract_bnd_mps(      Scalar * res_ptr      , std::array<long,3> res_dims,
                                                   const Scalar * const bnd_ptr, std::array<long,1> bnd_dims,
                                                   const Scalar * const mps_ptr, std::array<long,3> mps_dims){
//    auto t_con = tid::tic_token("contract_bnd_mps", tid::level::highest);
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,res_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(mps_ptr,mps_dims);
    auto bnd = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>(bnd_ptr,bnd_dims);
    if(mps_dims[1] != bnd_dims[0]) throw except::runtime_error("Dimension mismatch mps {} (idx 1) and bnd {} (idx 0)", mps_dims, bnd_dims);
    if(mps_dims != res_dims) throw except::runtime_error("Dimension mismatch mps {} and res {}", mps_dims, res_dims);
    auto &threads = tenx::threads::get();
    res.device(*threads->dev) = tenx::asDiagonal(bnd).contract(mps, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2});
}

template void tools::common::contraction::contract_bnd_mps(      cx64 *       res_ptr, std::array<long,3> res_dims,
                                                           const cx64 * const bnd_ptr, std::array<long,1> bnd_dims,
                                                           const cx64 * const mps_ptr, std::array<long,3> mps_dims);
template void tools::common::contraction::contract_bnd_mps(      fp64 *       res_ptr, std::array<long,3> res_dims,
                                                           const fp64 * const bnd_ptr, std::array<long,1> bnd_dims,
                                                           const fp64 * const mps_ptr, std::array<long,3> mps_dims);


template<typename Scalar>
void tools::common::contraction::contract_mps_mps(      Scalar * res_ptr       , std::array<long,3> res_dims,
                                                  const Scalar * const mpsL_ptr, std::array<long,3> mpsL_dims,
                                                  const Scalar * const mpsR_ptr, std::array<long,3> mpsR_dims){
//    auto t_con = tid::tic_token("contract_mps_mps", tid::level::highest);
    auto res  = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,res_dims);
    auto mpsL = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(mpsL_ptr, mpsL_dims);
    auto mpsR = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(mpsR_ptr, mpsR_dims);
    constexpr auto shuffle_idx  = std::array<long,4>{0, 2, 1, 3};
    constexpr auto contract_idx = tenx::idx({2}, {1});
    auto check_dims = std::array<long,3>{mpsL.dimension(0) * mpsR.dimension(0), mpsL.dimension(1), mpsR.dimension(2)};
    if(res_dims != check_dims) throw except::runtime_error("res dimension mismatch: dims {} | expected dims {}", res_dims, check_dims);
    if(mpsL.dimension(2) != mpsR.dimension(1)) throw except::runtime_error("Dimension mismatch mpsL {} (idx 2) and mpsR {} (idx 1)", mpsL.dimensions(), mpsR.dimensions());
    auto &threads = tenx::threads::get();
    if constexpr(std::is_same_v<Scalar, fp64>){
        auto tmp = Eigen::Tensor<Scalar,4>(mpsL_dims[0], mpsL_dims[1], mpsR_dims[0], mpsR_dims[2]);
        #if defined(DMRG_ENABLE_TBLIS)
        auto arch =  get_tblis_arch();
        const tblis::tblis_config_s *tblis_config = tblis::tblis_get_config(arch.data());
        contract_tblis(mpsL, mpsR, tmp, "abe", "ced", "abcd", tblis_config);
        res.device(*threads->dev)  = tmp.shuffle(shuffle_idx).reshape(res_dims);
        #else
        res.device(*threads->dev) = mpsL.contract(mpsR, contract_idx).shuffle(shuffle_idx).reshape(res_dims);
        #endif
    }else{
        res.device(*threads->dev) = mpsL.contract(mpsR, contract_idx).shuffle(shuffle_idx).reshape(res_dims);
    }
}


template void tools::common::contraction::contract_mps_mps(      cx64 * res_ptr       , std::array<long,3> res_dims,
                                                           const cx64 * const mpsL_ptr, std::array<long,3> mpsL_dims,
                                                           const cx64 * const mpsR_ptr, std::array<long,3> mpsR_dims);

template void tools::common::contraction::contract_mps_mps(      fp64 * res_ptr       , std::array<long,3> res_dims,
                                                           const fp64 * const mpsL_ptr, std::array<long,3> mpsL_dims,
                                                           const fp64 * const mpsR_ptr, std::array<long,3> mpsR_dims);



template<typename Scalar>
Scalar tools::common::contraction::contract_mps_mps_overlap(const Scalar * const mps1_ptr, std::array<long,3> mps1_dims,
                                                            const Scalar * const mps2_ptr, std::array<long,3> mps2_dims){
    auto t_con = tid::tic_token("contract_mps_mps_overlap", tid::level::highest);
    auto size1 = Eigen::internal::array_prod(mps1_dims);
    auto size2 = Eigen::internal::array_prod(mps2_dims);
    auto mps1 = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic,1>>(mps1_ptr, size1);
    auto mps2 = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic,1>>(mps2_ptr, size2);
    if(size1 != size2) throw except::runtime_error("Size mismatch mps1 {} and mps2 {}",size1, size2);
    return mps1.dot(mps2); // Calls gemv // TODO: Check that this works with the tests (used to be conjugate on mps2!)
}

template cx64   tools::common::contraction::contract_mps_mps_overlap(const cx64 * const mps1_ptr, std::array<long,3> mps1_dims,
                                                                     const cx64 * const mps2_ptr, std::array<long,3> mps2_dims);
template fp64   tools::common::contraction::contract_mps_mps_overlap(const fp64 * const mps1_ptr, std::array<long,3> mps1_dims,
                                                                     const fp64 * const mps2_ptr, std::array<long,3> mps2_dims);


template<typename Scalar>
void tools::common::contraction::contract_mps_mps_partial(Scalar *       res_ptr , std::array<long,2> res_dims,
                                                          const Scalar * const mps1_ptr, std::array<long,3> mps1_dims,
                                                          const Scalar * const mps2_ptr, std::array<long,3> mps2_dims,
                                                          std::array<long,2> idx){
    auto res  = Eigen::TensorMap<Eigen::Tensor<Scalar,2>>(res_ptr,res_dims);
    auto mps1 = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(mps1_ptr, mps1_dims);
    auto mps2 = Eigen::TensorMap<const Eigen::Tensor<Scalar,3>>(mps2_ptr, mps2_dims);
    auto idxs = tenx::idx(idx,idx);
    auto &threads = tenx::threads::get();
    res.device(*threads->dev) = mps1.conjugate().contract(mps2, idxs);
}

template void tools::common::contraction::contract_mps_mps_partial(      cx64 *       res_ptr , std::array<long,2> res_dims,
                                                                   const cx64 * const mps1_ptr, std::array<long,3> mps1_dims,
                                                                   const cx64 * const mps2_ptr, std::array<long,3> mps2_dims,
                                                                   std::array<long,2> idx);
template void tools::common::contraction::contract_mps_mps_partial(      fp64 *       res_ptr , std::array<long,2> res_dims,
                                                                   const fp64 * const mps1_ptr, std::array<long,3> mps1_dims,
                                                                   const fp64 * const mps2_ptr, std::array<long,3> mps2_dims,
                                                                   std::array<long,2> idx);


template<typename Scalar>
void tools::common::contraction::contract_env_mps_mpo(      Scalar *      res_ptr, std::array<long, 2> res_dims,
                                                      const Scalar *const env_ptr, std::array<long, 2> env_dims,
                                                      const Scalar *const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar *const mpo_ptr, std::array<long, 2> mpo_dims) {
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar, 2>>(res_ptr, res_dims);
    auto env = Eigen::TensorMap<const Eigen::Tensor<Scalar, 2>>(env_ptr, env_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<Scalar, 2>>(mpo_ptr, mpo_dims);
    auto &threads = tenx::threads::get();
    res.device(*threads->dev) = env.contract(mps,             tenx::idx({0}, {1}))
                                  .contract(mpo,             tenx::idx({1}, {0}))
                                  .contract(mps.conjugate(), tenx::idx({0, 2}, {1, 0}));
}

template void tools::common::contraction::contract_env_mps_mpo(      cx64 *       res_ptr , std::array<long,2> res_dims,
                                                               const cx64 * const env_ptr , std::array<long,2> env_dims,
                                                               const cx64 * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cx64 * const mpo_ptr , std::array<long,2> mpo_dims);

template<typename Scalar>
void tools::common::contraction::contract_env_mps_mpo(      Scalar *       res_ptr, std::array<long, 3> res_dims,
                                                      const Scalar * const env_ptr, std::array<long, 3> env_dims,
                                                      const Scalar * const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar * const mpo_ptr, std::array<long, 4> mpo_dims) {
    auto res                           = Eigen::TensorMap<Eigen::Tensor<Scalar, 3>>(res_ptr, res_dims);
    auto env                           = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(env_ptr, env_dims);
    auto mps                           = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo                           = Eigen::TensorMap<const Eigen::Tensor<Scalar, 4>>(mpo_ptr, mpo_dims);
    auto &threads = tenx::threads::get();
    res.device(*threads->dev) = env.contract(mps,             tenx::idx({0}, {1}))
                                  .contract(mpo,             tenx::idx({1, 2}, {0, 2}))
                                  .contract(mps.conjugate(), tenx::idx({0, 3}, {1, 0}))
                                  .shuffle(                  tenx::array3{0, 2, 1});
}

template void tools::common::contraction::contract_env_mps_mpo(      cx64 *       res_ptr , std::array<long,3> res_dims,
                                                               const cx64 * const env_ptr , std::array<long,3> env_dims,
                                                               const cx64 * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cx64 * const mpo_ptr , std::array<long,4> mpo_dims);

template<typename Scalar>
void tools::common::contraction::contract_mps_mpo_env(      Scalar *       res_ptr, std::array<long, 2> res_dims,
                                                      const Scalar * const env_ptr, std::array<long, 2> env_dims,
                                                      const Scalar * const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar * const mpo_ptr, std::array<long, 2> mpo_dims) {
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar, 2>>(res_ptr, res_dims);
    auto env = Eigen::TensorMap<const Eigen::Tensor<Scalar, 2>>(env_ptr, env_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<Scalar, 2>>(mpo_ptr, mpo_dims);
    auto &threads = tenx::threads::get();
    res.device(*threads->dev) = env.contract(mps,             tenx::idx({0}, {2}))
                                  .contract(mpo,             tenx::idx({1}, {0}))
                                  .contract(mps.conjugate(), tenx::idx({0, 2}, {2, 0}));
}
template void tools::common::contraction::contract_mps_mpo_env(      cx64 *       res_ptr , std::array<long,2> res_dims,
                                                               const cx64 * const env_ptr , std::array<long,2> env_dims,
                                                               const cx64 * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cx64 * const mpo_ptr , std::array<long,2> mpo_dims);
template<typename Scalar>
void tools::common::contraction::contract_mps_mpo_env(Scalar       *      res_ptr, std::array<long, 3> res_dims,
                                                      const Scalar * const env_ptr, std::array<long, 3> env_dims,
                                                      const Scalar * const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar * const mpo_ptr, std::array<long, 4> mpo_dims) {
    auto res                           = Eigen::TensorMap<Eigen::Tensor<Scalar, 3>>(res_ptr, res_dims);
    auto env                           = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(env_ptr, env_dims);
    auto mps                           = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo                           = Eigen::TensorMap<const Eigen::Tensor<Scalar, 4>>(mpo_ptr, mpo_dims);
    auto &threads = tenx::threads::get();
    res.device(*threads->dev) = env.contract(mps,             tenx::idx({0}, {2}))
                                  .contract(mpo,             tenx::idx({1, 2}, {1, 2}))
                                  .contract(mps.conjugate(), tenx::idx({0, 3}, {2, 0}))
                                  .shuffle(                  tenx::array3{0, 2, 1});
}

template void tools::common::contraction::contract_mps_mpo_env(      cx64 *       res_ptr , std::array<long,3> res_dims,
                                                               const cx64 * const env_ptr , std::array<long,3> env_dims,
                                                               const cx64 * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cx64 * const mpo_ptr , std::array<long,4> mpo_dims);
/* clang-format on */
