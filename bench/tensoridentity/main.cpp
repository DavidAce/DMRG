#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"
#include "env/environment.h"
#include <fmt/core.h>
#include <string_view>
#include <unsupported/Eigen/CXX11/Tensor>

template<typename Scalar = double>
Eigen::Tensor<Scalar, 2> identity1(const Eigen::Index &dim) {
    using M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    M id    = M::Identity(dim, dim);
    return Eigen::TensorMap<const Eigen::Tensor<double, 2>>(id.data(), dim, dim).cast<Scalar>();
}

template<typename Scalar = double>
Eigen::Tensor<Scalar, 2> identity2(const Eigen::Index &dim) {
    Eigen::Tensor<double, 1> tensor(dim);
    tensor.setConstant(1);
    return tensor.inflate(std::array<long, 1>{tensor.size() + 1}).reshape(std::array<long, 2>{tensor.size(), tensor.size()}).template cast<Scalar>();
}

int main() {
    fmt::print("Compiler flags {}", env::build::compiler_flags);
    using real = double;
    using cplx = std::complex<double>;
    Eigen::Tensor<real, 2>        idd;
    Eigen::Tensor<cplx, 2>        idc;
    constexpr std::array<long, 5> dims{200l, 500l, 5000l};

    ankerl::nanobench::Bench().run("real identity1", [&] {
        for(const auto &dim : dims) { idd = identity1<real>(dim); }
        ankerl::nanobench::doNotOptimizeAway(idd);
    });
    ankerl::nanobench::Bench().run("cplx identity1", [&] {
        for(const auto &dim : dims) { idc = identity1<cplx>(dim); }
        ankerl::nanobench::doNotOptimizeAway(idc);
    });
    ankerl::nanobench::Bench().run("real identity2", [&] {
        for(const auto &dim : dims) { idd = identity2<real>(dim); }
        ankerl::nanobench::doNotOptimizeAway(idd);
    });
    ankerl::nanobench::Bench().run("cplx identity2", [&] {
        for(const auto &dim : dims) { idc = identity2<cplx>(dim); }
        ankerl::nanobench::doNotOptimizeAway(idc);
    });
}
