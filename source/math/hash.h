#pragma once
namespace hash{
    template<typename T>
    std::size_t hash_buffer(const T * v, unsigned long size, std::size_t seed = 0);
}
