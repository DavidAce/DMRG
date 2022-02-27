#pragma once
#include "tools/common/log.h"
#include <optional>
#include <string>
#include <vector>
namespace tools::finite::opt {
    class opt_mps;
}

namespace tools::finite::opt::reports {
    struct bfgs_entry {
        std::string description;
        long        size, space;
        double      energy, variance, overlap, norm;
        double      delta_f, max_grad_norm;
        size_t      iter, counter;
        double      time;
    };
    struct time_entry {
        double vH2v, vHv, vH2, vH, bfgs, step;
    };
    struct subs_entry {
        long   nev;
        double max_olap, min_olap, eps, eig_time, ham_time, lu_time;
        long   iter, mv, pc;
    };

    struct eigs_entry {
        std::string               description;
        std::string               ritz;
        long                      size, idx, nev, ncv;
        double                    energy, eigval, variance, overlap, norm, tol, resid, grad;
        size_t                    iter, mv, pc;
        double                    time;
        spdlog::level::level_enum level = spdlog::level::debug;
    };

    inline std::vector<bfgs_entry> bfgs_log;
    inline std::vector<time_entry> time_log;
    inline std::vector<subs_entry> subs_log;
    inline std::vector<eigs_entry> eigs_log;

    extern void print_bfgs_report();
    extern void print_time_report();
    extern void print_subs_report();
    extern void print_eigs_report(std::optional<size_t> max_entries = std::nullopt);
    extern void bfgs_add_entry(const std::string &description, long size, long space, double energy, double variance, double overlap, double norm,
                               double delta_f, double grad_norm, size_t iter, size_t counter, double time);
    extern void bfgs_add_entry(std::string_view mode, std::string_view tag, const opt_mps &mps, std::optional<long> space = std::nullopt);
    extern void time_add_entry();
    extern void subs_add_entry(long nev, double max_olap, double min_olap, double eps, double eig_time, double ham_time, double lu_time, long iter, long mv,
                               long pc);
    extern void eigs_add_entry(const opt_mps &mps, spdlog::level::level_enum level = spdlog::level::debug);
}