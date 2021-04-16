#pragma once
#include <fmt/core.h>
namespace svd {
    struct settings {
        std::optional<double> threshold   = std::nullopt;
        std::optional<long>   switchsize  = std::nullopt;
        std::optional<size_t> loglevel    = std::nullopt;
        std::optional<bool>   use_bdc     = std::nullopt;
        std::optional<bool>   use_lapacke = std::nullopt;
        std::optional<bool>   profile     = std::nullopt;
        /* clang-format off */
        std::string to_string(){
            std::string msg = "svd";
            if(threshold   ) msg.append(fmt::format(" | threshold {:.2e}", threshold.value()));
            if(switchsize  ) msg.append(fmt::format(" | switchsize {}", switchsize.value()));
            if(loglevel    ) msg.append(fmt::format(" | loglevel {}", loglevel.value()));
            if(use_bdc     ) msg.append(fmt::format(" | use_bdc {}", use_bdc.value()));
            if(use_lapacke ) msg.append(fmt::format(" | use_lapacke {}", use_lapacke.value()));
            if(profile     ) msg.append(fmt::format(" | profile {}", profile.value()));
            return msg;
        }
        /* clang-format on */
    };
}