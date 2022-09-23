#pragma once
#include <fmt/core.h>
#include <optional>
#include <string>

namespace svd {
    enum class lib { eigen, lapacke, rsvd };
    inline std::string_view enum2sv(svd::lib lib) {
        if(lib == svd::lib::eigen) return "eigen";
        if(lib == svd::lib::lapacke) return "lapacke";
        if(lib == svd::lib::rsvd) return "rsvd";
        throw std::logic_error("Could not match svd::lib");
    }

    struct config {
        std::optional<long>     rank_max       = std::nullopt;
        std::optional<long>     rank_min       = std::nullopt; /*!< Keep this many singular values even if they are smaller than truncation_lim */
        std::optional<double>   truncation_lim = std::nullopt;
        std::optional<size_t>   switchsize_bdc = std::nullopt;
        std::optional<size_t>   loglevel       = std::nullopt;
        std::optional<bool>     use_bdc        = std::nullopt;
        std::optional<svd::lib> svd_lib        = std::nullopt;
        std::optional<bool>     save_fail      = std::nullopt;
        std::optional<bool>     benchmark      = std::nullopt;
        /* clang-format off */
        std::string to_string(){
            std::string msg;
            if(rank_max        ) msg.append(fmt::format(" | rank_max {}", rank_max.value()));
            if(rank_min        ) msg.append(fmt::format(" | rank_min {}", rank_min.value()));
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
        config() = default;
        explicit config(long rank_max_);
        explicit config(double truncation_lim_);
        explicit config(long rank_max_, double truncation_lim_);
        explicit config(std::optional<long> rank_max_);
        explicit config(std::optional<double> truncation_lim_);
        explicit config(std::optional<long> rank_max_, std::optional<double> truncation_lim_);
    };
}