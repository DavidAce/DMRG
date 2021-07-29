#pragma once

#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#if defined(FMT_HEADER_ONLY)
    #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
#endif