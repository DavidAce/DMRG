#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"
#include "tid/tid.h"
#include <fmt/core.h>
#include <string_view>

void scope() { auto t_tid = tid::tic_scope("scope"); }
void token() { auto t_tid = tid::tic_token("token"); }
int  main() {
    tid::set_level(tid::level::normal);
    auto t_total = tid::tic_scope("total", tid::level::higher);
    auto bench   = ankerl::nanobench::Bench().title("tid").relative(true).epochs(10000).minEpochIterations(10000);
    ankerl::nanobench::Bench().run("scope", [&] { scope(); });
    ankerl::nanobench::Bench().run("token", [&] { token(); });
    ankerl::nanobench::Bench().run("many scopes", [&] {
        auto t_scope1 = tid::tic_scope("scope1", tid::level::higher);
        auto t_scope2 = tid::tic_scope("scope2", tid::level::higher);
        auto t_scope3 = tid::tic_scope("scope3", tid::level::higher);
        auto t_scope4 = tid::tic_scope("scope4", tid::level::higher);
        ankerl::nanobench::doNotOptimizeAway(t_scope1);
        ankerl::nanobench::doNotOptimizeAway(t_scope2);
        ankerl::nanobench::doNotOptimizeAway(t_scope3);
        ankerl::nanobench::doNotOptimizeAway(t_scope4);
    });
    for(const auto &t : tid::get_tree("total")) fmt::print("{}\n", t.str());

}
