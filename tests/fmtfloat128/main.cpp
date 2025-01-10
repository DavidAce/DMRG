
#include "io/fmt_f128_t.h"
#include "math/f128.h"
#include <fmt/format.h>
#include <cassert>
#include <fmt/std.h>


int main() {
    f128_t v1 = f128_t("-3721.1244523527272214135626241245632");
    fp128  v2 = f128_t("3721.1244523527272214135626241245632");
    fp128  v3 = 3;
    fp128  v4 = 3.14;
    auto   v5 = cx128(3.14, 2.3);
    cx128  v6 = 0;
    cx128  v7 = v1.value() * v6;

    fmt::print("v1: (:.31f)     : [{:.31f}]\n", v1);
    fmt::print("v1: (:>+50.31f) : [{:>+50.31f}]\n", v1);
    fmt::print("v1: (:<+50.31f) : [{:<+50.31f}]\n", v1);
    auto str1a =  fmt::format("{:.31f}", v1);
    auto str1b =  fmt::format("{:>+50.31f}", v1);
    auto str1c =  fmt::format("{:<+50.31f}", v1);

    fmt::print("str1a : [{}] (size: {})\n", str1a, str1a.size());
    assert(str1a == "-3721.1244523527272214135626241245632");
    assert(str1b == "             -3721.1244523527272214135626241245632");
    assert(str1c == "-3721.1244523527272214135626241245632             ");



    fmt::print("v2: (:.31f)     : [{:.31f}]\n",    f128_t(v2));
    fmt::print("v2: (:+>50.31f) : [{:>+50.31f}]\n", f128_t(v2));
    fmt::print("v2: (:+<50.31f) : [{:<+50.31f}]\n", f128_t(v2));
    auto str2a =  fmt::format("{:.31f}", f128_t(v2));
    auto str2b =  fmt::format("{:>+50.31f}", f128_t(v2));
    auto str2c =  fmt::format("{:<+50.31f}", f128_t(v2));

    assert(str2a == "3721.1244523527272214135626241245632");
    assert(str2b == "             +3721.1244523527272214135626241245632");
    assert(str2c == "+3721.1244523527272214135626241245632             ");
    fmt::print("v3: {:>20.7f}\n", f128_t(v3));
    fmt::print("v4: {:<20.7f}\n", f128_t(v4));
    fmt::print("v4.string(5,8,f,<): {}\n", f128_t(v4).string(5, 8, 'f', "<"));


    fmt::print("v5: {:>.5f}{:>+.5f}i\n", c128_t(v5).real(), c128_t(v5).imag());
    fmt::print("v6: {:>.5f}{:>+.5f}i\n", c128_t(v6).real(), c128_t(v6).imag());
    fmt::print("v7: {:>.5f}{:>+.5f}i\n", c128_t(v7).real(), c128_t(v7).imag());
    return 0;
}
