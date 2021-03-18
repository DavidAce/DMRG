#pragma once

namespace svd{
    struct settings{
        std::optional<double>   threshold   = std::nullopt;
        std::optional<long>     switchsize  = std::nullopt;
        std::optional<size_t>   loglevel    = std::nullopt;
        std::optional<bool>     use_bdc     = std::nullopt;
        std::optional<bool>     use_lapacke = std::nullopt;
        std::optional<bool>     profile     = std::nullopt;
    };
}