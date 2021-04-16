#pragma once
#include <optional>
#include <string>
#include <vector>
namespace tools::finite::opt {
    class opt_mps;
}

namespace tools::finite::opt::internal::reports {
    struct bfgs_entry {
        std::string description;
        long        size, space;
        double      energy, variance, overlap, norm;
        double      delta_f, grad_max_norm;
        size_t      iter, counter;
        double      time;
    };
    struct time_entry {
        double vH2v, vHv, vH2, vH, bfgs, step;
    };
    struct eigs_entry {
        long   nev;
        double max_olap, min_olap, eps, eig_time, ham_time, lu_time;
        size_t steps;
    };

    struct krylov_entry {
        std::string description,ritz;
        long        size, nev, ncv;
        double      energy, eigval, variance, overlap, norm, tol;
        size_t      iter, counter;
        double      time;
    };

    inline std::vector<bfgs_entry>   bfgs_log;
    inline std::vector<time_entry>   time_log;
    inline std::vector<eigs_entry>   eigs_log;
    inline std::vector<krylov_entry> krylov_log;

    extern void print_bfgs_report();
    extern void print_time_report();
    extern void print_eigs_report();
    extern void print_krylov_report(std::optional<size_t> max_entries = std::nullopt);
    extern void bfgs_add_entry(const std::string &description, long size, long space, double energy, double variance, double overlap, double norm,
                               double delta_f, double grad_norm, size_t iter, size_t counter, double time);
    extern void bfgs_add_entry(const std::string &mode, const std::string &tag, const opt_mps &mps, std::optional<long> space = std::nullopt);
    extern void time_add_dir_entry();
    extern void time_add_sub_entry();
    extern void eigs_add_entry(long nev, double max_olap, double min_olap, double eps, double eig_time, double ham_time, double lu_time, size_t steps);
    extern void krylov_add_entry(const opt_mps &mps);
}