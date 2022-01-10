#pragma once
#include <fmt/core.h>
#include <optional>
#include <string>

enum class SVDLib { eigen, lapacke, rsvd };

namespace svd {
    inline std::string_view enum2sv(SVDLib lib) {
        if(lib == SVDLib::eigen) return "eigen";
        if(lib == SVDLib::lapacke) return "lapacke";
        if(lib == SVDLib::rsvd) return "rsvd";
        throw std::logic_error("Could not match SVDLib");
    }

    struct settings {
        std::optional<double> threshold      = std::nullopt;
        std::optional<size_t> switchsize_bdc = std::nullopt;
        std::optional<size_t> loglevel       = std::nullopt;
        std::optional<bool>   use_bdc        = std::nullopt;
        std::optional<SVDLib> svd_lib        = std::nullopt;
        std::optional<bool>   save_fail      = std::nullopt;
        std::optional<bool>   benchmark      = std::nullopt;
        /* clang-format off */
        std::string to_string(){
            std::string msg;
            if(threshold       ) msg.append(fmt::format(" | threshold {:.2e}", threshold.value()));
            if(switchsize_bdc  ) msg.append(fmt::format(" | switchsize bdc {}", switchsize_bdc.value()));
            if(loglevel        ) msg.append(fmt::format(" | loglevel {}", loglevel.value()));
            if(use_bdc         ) msg.append(fmt::format(" | use_bdc {}", use_bdc.value()));
            if(svd_lib         ) msg.append(fmt::format(" | svd_lib {}", enum2sv(svd_lib.value())));
            if(save_fail       ) msg.append(fmt::format(" | save_fail {}",save_fail.value()));
            if(benchmark       ) msg.append(fmt::format(" | benchmark {}",benchmark.value()));
            return msg.empty() ? msg : "svd settings" + msg;
        }
        /* clang-format on */
    };
}