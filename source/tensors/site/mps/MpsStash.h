#pragma once

template<typename T>
class stash {
    public:
    T      data;
    double error = 0;
    size_t pos_dst{};
};