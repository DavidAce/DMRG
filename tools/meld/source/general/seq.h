#pragma once
#include <algorithm>
#include <optional>
namespace seq {
    template<typename T>
    [[nodiscard]] std::optional<std::size_t> has(const T &sequence, const typename T::value_type &element) {
        if(std::empty(sequence)) return std::nullopt;
        auto it = std::find(std::begin(sequence), std::end(sequence), element);
        if(it == std::end(sequence))
            return std::nullopt;
        else
            return static_cast<size_t>(std::distance(std::begin(sequence), it));
    }

}