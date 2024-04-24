#include "human.h"
#include <array>
#include <fmt/compile.h>
#include <fmt/core.h>
#include <string_view>

namespace tools {
    std::string fmtBytes(bool on, size_t bytes, size_t base, size_t decimals) {
        if(on) {
            if(base != 1000 and base != 1024) throw std::runtime_error("Base must be either 1000 or 1024");
            constexpr uint64_t                             orders        = 6;
            constexpr std::array<std::string_view, orders> suffixes_1024 = {"B  ", "KiB", "MiB", "GiB", "TiB", "PiB"};
            constexpr std::array<std::string_view, orders> suffixes_1000 = {"B ", "kB", "MB", "GB", "TB", "PB"};
            auto                                           dblBase       = static_cast<double>(base);
            auto                                           dblByte       = static_cast<double>(bytes);

            uint64_t i;
            for(i = 0; i < orders and dblByte >= dblBase; i++, bytes /= base) dblByte = static_cast<double>(bytes) / dblBase;
            if(base == 1024)
                return fmt::format(FMT_COMPILE("{:.{}f} {}"), dblByte, decimals, suffixes_1024[i]);
            else
                return fmt::format(FMT_COMPILE("{:.{}f} {}"), dblByte, decimals, suffixes_1000[i]);
        } else
            return fmt::format(FMT_COMPILE("{}"), bytes);
    }
}