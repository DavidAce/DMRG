#pragma once
#include <fmt/core.h>
#include <optional>
#include <string>

namespace svd {
    enum class Lib { eigen, lapacke, rsvd };
    inline std::string_view enum2sv(svd::Lib lib) {
        if(lib == svd::Lib::eigen) return "eigen";
        if(lib == svd::Lib::lapacke) return "lapacke";
        if(lib == svd::Lib::rsvd) return "rsvd";
        throw std::logic_error("Could not match SVDLib");
    }

    struct settings {
        std::optional<long>     rank_max       = std::nullopt;
        std::optional<double>   truncation_lim = std::nullopt;
        std::optional<size_t>   switchsize_bdc = std::nullopt;
        std::optional<size_t>   loglevel       = std::nullopt;
        std::optional<bool>     use_bdc        = std::nullopt;
        std::optional<svd::Lib> svd_lib        = std::nullopt;
        std::optional<bool>     save_fail      = std::nullopt;
        std::optional<bool>     benchmark      = std::nullopt;
        /* clang-format off */
        std::string to_string(){
            std::string msg;
            if(truncation_lim  ) msg.append(fmt::format(" | truncation_lim {:.2e}", truncation_lim.value()));
            if(switchsize_bdc  ) msg.append(fmt::format(" | switchsize bdc {}", switchsize_bdc.value()));
            if(loglevel        ) msg.append(fmt::format(" | loglevel {}", loglevel.value()));
            if(use_bdc         ) msg.append(fmt::format(" | use_bdc {}", use_bdc.value()));
            if(svd_lib         ) msg.append(fmt::format(" | svd_lib {}", enum2sv(svd_lib.value())));
            if(save_fail       ) msg.append(fmt::format(" | save_fail {}",save_fail.value()));
            if(benchmark       ) msg.append(fmt::format(" | benchmark {}",benchmark.value()));
            return msg.empty() ? msg : "svd settings" + msg;
        }
        /* clang-format on */
        settings() = default;
        explicit settings(long rank_max_);
        explicit settings(double truncation_lim_);
        explicit settings(long rank_max_, double truncation_lim_);
    };
}