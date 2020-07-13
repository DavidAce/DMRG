#pragma once

namespace tools::common::moments{
    template<typename mps_type, typename mpo_type, typename env_type>
    extern double first(const mps_type &mps,const mpo_type &mpo, const env_type &envL, const env_type &envR);
    template<typename mps_type, typename mpo_type, typename env_type>
    extern double second(const mps_type &mps,const mpo_type &mpo, const env_type &envL, const env_type &envR);
}