#define ANKERL_NANOBENCH_IMPLEMENT
#include "config/settings.h"
#include "config/threading.h"
#include "env/environment.h"
#include "math/tenx.h"
#include "nanobench.h"
#include "tools/common/log.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <string_view>
#include <tools/common/contraction.h>
#include <unsupported/Eigen/CXX11/Tensor>

#if defined(DMRG_ENABLE_TBLIS)
    #include <tblis/tblis.h>
    #include <tblis/util/thread.h>
    #include <tci/tci_config.h>
#endif

using namespace tools::common::contraction;
template<typename res_type, typename mps_type, typename mpo_type, typename env_type>
void matrix_vector_product_eigen(TensorWrite<res_type> &res_, const TensorRead<mps_type> &mps_, const TensorRead<mpo_type> &mpo_,
                                 const TensorRead<env_type> &enL_, const TensorRead<env_type> &enR_) {
    using T       = typename mps_type::Scalar;
    auto  res     = tenx::asEval(res_);
    auto  mps     = tenx::asEval(mps_);
    auto  mpo     = tenx::asEval(mpo_);
    auto  enL     = tenx::asEval(enL_);
    auto  enR     = tenx::asEval(enR_);
    auto &threads = tenx::threads::get();

    res.device(*threads->dev) =
        mps.contract(enL, tenx::idx({1}, {0})).contract(mpo, tenx::idx({3, 0}, {0, 2})).contract(enR, tenx::idx({0, 2}, {0, 2})).shuffle(tenx::array3{1, 0, 2});
}

template<typename res_type, typename mps_type, typename mpo_type, typename env_type>
void matrix_vector_product_eigen2(TensorWrite<res_type> &res_, const TensorRead<mps_type> &mps_, const TensorRead<mpo_type> &mpo_,
                                  const TensorRead<env_type> &enL_, const TensorRead<env_type> &enR_) {
    using T       = typename mps_type::Scalar;
    auto  res     = tenx::asEval(res_);
    auto  mps     = tenx::asEval(mps_);
    auto  mpo     = tenx::asEval(mpo_);
    auto  enL     = tenx::asEval(enL_);
    auto  enR     = tenx::asEval(enR_);
    auto &threads = tenx::threads::get();

    res.device(*threads->dev) = mps.contract(enL, tenx::idx({0}, {0})).contract(mpo, tenx::idx({2, 0}, {0, 2})).contract(enR, tenx::idx({0, 2}, {0, 1}));
}

#if defined(DMRG_ENABLE_TBLIS)
template<typename ea_type, typename eb_type, typename ec_type>
void contract_tblis(const TensorRead<ea_type> &ea, const TensorRead<eb_type> &eb, TensorWrite<ec_type> &ec, const tblis::label_vector &la,
                    const tblis::label_vector &lb, const tblis::label_vector &lc, std::string_view arch) {
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

    const tblis::tblis_config_s *tblis_config = tblis::tblis_get_config(arch.data());
    #if defined(TCI_USE_OPENMP_THREADS) && defined(_OPENMP)
    tblis_set_num_threads(static_cast<unsigned int>(omp_get_max_threads()));
    #endif
    tblis_tensor_mult(nullptr, tblis_config, &A_s, la.c_str(), &B_s, lb.c_str(), &C_s, lc.c_str());
}


template<typename res_type, typename mps_type, typename mpo_type, typename env_type>
void matrix_vector_product_tblis(TensorWrite<res_type> &res_, const TensorRead<mps_type> &mps_, const TensorRead<mpo_type> &mpo_,
                                 const TensorRead<env_type> &enL_, const TensorRead<env_type> &enR_, std::string arch) {
    using T  = typename mps_type::Scalar;
    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto  res     = tenx::asEval(res_);
    auto  mps     = tenx::asEval(mps_);
    auto  mpo     = tenx::asEval(mpo_);
    auto  enL     = tenx::asEval(enL_);
    auto  enR     = tenx::asEval(enR_);
    auto &threads = tenx::threads::get();

    long d       = mps.dimension(0);
    long mL      = mpo.dimension(0);
    long mR      = mpo.dimension(1);
    long chiL    = mps.dimension(1);
    long chiR    = mps.dimension(2);
    auto mps_mat = Eigen::Map<const MT>(mps.data(), d * chiL, chiR);
    auto enR_mat = Eigen::Map<const MT>(enR.data(), chiR, chiR * mR);


    Eigen::Tensor<T, 4> mpsenvL(mps.dimension(0), mps.dimension(2), enL.dimension(1), enL.dimension(2));
    Eigen::Tensor<T, 4> mpsenvLmpo(mps.dimension(2), enL.dimension(1), mpo.dimension(1), mpo.dimension(3));
    contract_tblis(mps, enL, mpsenvL, "afb", "fcd", "abcd", arch);
    contract_tblis(mpsenvL, mpo, mpsenvLmpo, "qijr", "rkql", "ijkl", arch);
    contract_tblis(mpsenvLmpo, enR, res, "qjri", "qkr", "ijk", arch);
}
#endif

int main() {
    tools::log = tools::Logger::setLogger("matvec", 2, true);
    fmt::print("Compiler flags {}", env::build::compiler_flags);
    settings::threading::num_threads = 4;
    settings::configure_threads();
    using real = double;
    using cplx = std::complex<double>;
    long d     = 8;
    long m     = 12;
    long chiL  = 128;
    long chiR  = 128;

    tools::log->info("__builtin_cpu_supports(avx512f)    : {}", __builtin_cpu_supports("avx512f"));
    tools::log->info("__builtin_cpu_supports(x86-64-v4)  : {}", __builtin_cpu_supports("x86-64-v4"));
    tools::log->info("__builtin_cpu_supports(x86-64-v3)  : {}", __builtin_cpu_supports("x86-64-v3"));
    tools::log->info("__builtin_cpu_supports(x86-64-v2)  : {}", __builtin_cpu_supports("x86-64-v2"));
    tools::log->info("__builtin_cpu_supports(x86-64)     : {}", __builtin_cpu_supports("x86-64"));
    tools::log->info("__builtin_cpu_is(znver1)            : {}", __builtin_cpu_is("znver1"));
    tools::log->info("__builtin_cpu_is(znver2)            : {}", __builtin_cpu_is("znver2"));
    tools::log->info("__builtin_cpu_is(znver3)            : {}", __builtin_cpu_is("znver3"));
    tools::log->info("__builtin_cpu_is(znver4)            : {}", __builtin_cpu_is("znver4"));
    // tools::log->info("__builtin_cpu_is(znver5)            : {}", __builtin_cpu_is("znver5"));

    {
        auto res = Eigen::Tensor<real, 3>(d, chiL, chiR);
        auto mps = Eigen::Tensor<real, 3>(d, chiL, chiR);
        auto mpo = Eigen::Tensor<real, 4>(m, m, d, d);
        auto enL = Eigen::Tensor<real, 3>(chiL, chiL, m);
        auto enR = Eigen::Tensor<real, 3>(chiR, chiR, m);
        tools::log->info("mps: size {} | {}", mps.size(), mps.dimensions());

        mps.setRandom();
        mpo.setRandom();
        enL.setRandom();
        enR.setRandom();
        // ankerl::nanobench::Bench().warmup(5).epochs(40).run("tblis haswell", [&] {
        //     tools::common::contraction::matrix_vector_product(res, mps, mpo, enL, enR);
        //     ankerl::nanobench::doNotOptimizeAway(res);
        //     ankerl::nanobench::doNotOptimizeAway(mps);
        // });
#if defined(DMRG_ENABLE_TBLIS)
        ankerl::nanobench::Bench().warmup(5).epochs(40).run("tblis haswell", [&] {
            matrix_vector_product_tblis(res, mps, mpo, enL, enR, "haswell");
            ankerl::nanobench::doNotOptimizeAway(res);
            ankerl::nanobench::doNotOptimizeAway(mps);
        });
        ankerl::nanobench::Bench().warmup(5).epochs(40).run("tblis zen", [&] {
            matrix_vector_product_tblis(res, mps, mpo, enL, enR, "zen");
            ankerl::nanobench::doNotOptimizeAway(res);
            ankerl::nanobench::doNotOptimizeAway(mps);
        });
        // ankerl::nanobench::Bench().warmup(5).epochs(40).run("tblis skx1", [&] {
        //     matrix_vector_product_tblis(res, mps, mpo, enL, enR, "intel");
        //     ankerl::nanobench::doNotOptimizeAway(res);
        //     ankerl::nanobench::doNotOptimizeAway(mps);
        // });
        // ankerl::nanobench::Bench().warmup(5).epochs(40).run("tblis skx2", [&] {
        //     matrix_vector_product_tblis(res, mps, mpo, enL, enR, "skx2");
        //     ankerl::nanobench::doNotOptimizeAway(res);
        //     ankerl::nanobench::doNotOptimizeAway(mps);
        // });
#endif
    }

    {
        auto res = Eigen::Tensor<real, 3>(chiL, d, chiR);
        auto mps = Eigen::Tensor<real, 3>(chiL, d, chiR);
        auto mpo = Eigen::Tensor<real, 4>(m, m, d, d);
        auto enL = Eigen::Tensor<real, 3>(chiL, m, chiL);
        auto enR = Eigen::Tensor<real, 3>(chiR, m, chiR);

        mps.setRandom();
        mpo.setRandom();
        enL.setRandom();
        enR.setRandom();

        ankerl::nanobench::Bench().warmup(5).epochs(40).run("matvec eigen2", [&] {
            matrix_vector_product_eigen2(res, mps, mpo, enL, enR);
            ankerl::nanobench::doNotOptimizeAway(res);
            ankerl::nanobench::doNotOptimizeAway(mps);
        });
    }
}
