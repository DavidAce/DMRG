#include "algorithms/AlgorithmStatus.h"
#include "config/debug.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "math/tenx/threads.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <string>

tools::finite::opt::opt_mps tools::finite::opt::get_opt_initial_mps(const TensorsFinite &tensors) {
    auto    t_init = tid::tic_scope("initial_mps", tid::level::higher);
    opt_mps initial_mps("current mps", tensors.get_multisite_mps(), tensors.active_sites,
                        tools::finite::measure::energy_shift(tensors),              // Shifted energy for full system
                        tools::finite::measure::energy_minus_energy_shift(tensors), // Eigval
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap to initial state (self) is obviously 1
                        tensors.get_length());
    initial_mps.validate_initial_mps();
    return initial_mps;
}

tools::finite::opt::opt_mps tools::finite::opt::find_excited_state(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                                   OptMeta &meta) {
    auto t_opt        = tid::tic_scope("opt");
    meta.problem_size = tensors.active_problem_size();
    if(meta.optSolver == OptSolver::EIG and meta.problem_size > settings::precision::eig_max_size) {
        tools::log->debug("Problem size {} > eig_max_size ({}): Switching solver {} -> {}", meta.problem_size, settings::precision::eig_max_size,
                          enum2sv(meta.optSolver), enum2sv(OptSolver::EIGS));
        meta.optSolver = OptSolver::EIGS;
    }
    tools::log->trace("Starting optimization: func [{}] | algo [{}] | solver [{}] | type [{}] | ritz [{}] | position [{}] | sites {} | shape {} = {}",
                      enum2sv(meta.optCost), enum2sv(meta.optAlgo), enum2sv(meta.optSolver), enum2sv(meta.optType), enum2sv(meta.optRitz), status.position,
                      tensors.active_sites, tensors.active_problem_dims(), tensors.active_problem_size());

    using namespace opt::internal;

    if(initial_mps.get_sites() != tensors.active_sites)
        throw except::runtime_error("mismatch in active sites: initial_mps {} | active {}", initial_mps.get_sites(), tensors.active_sites);
    opt_mps result;
    switch(meta.optCost) {
        case OptCost::OVERLAP: result = internal::optimize_overlap(tensors, initial_mps, meta); break;
        case OptCost::ENERGY: result = internal::optimize_energy_eigs(tensors, initial_mps, status, meta); break;
        case OptCost::VARIANCE: {
            switch(meta.optAlgo) {
                case OptAlgo::DIRECT: {
                    result = internal::optimize_variance_eigs(tensors, initial_mps, status, meta);
                    break;
                }
                case OptAlgo::DIRECTX2: {
                    auto metax2     = meta;
                    metax2.optCost  = OptCost::ENERGY;
                    metax2.eigs_tol = std::max(1e-12, meta.eigs_tol.value_or(1e-12));
                    metax2.eigs_ncv = std::max(32, meta.eigs_ncv.value_or(32));
                    // metax2.primme_method = "PRIMME_DEFAULT_MIN_MATVECS";
                    result = internal::optimize_energy_eigs(tensors, initial_mps, status, metax2);
                    tools::finite::opt::reports::print_eigs_report();
                    if(meta.optRitz == OptRitz::SM) {
                        // We accept this result if the residual norm is small, and either
                        // - The energy is smaller in absolute
                        // - The variance is lower
                        // Otherwise we call optimize_variance_eigs on the result
                        // auto sigma    = std::sqrt(initial_mps.get_variance());
                        // bool in2sigma = std::abs(result.get_energy() - initial_mps.get_energy()) <= 2 * sigma;
                        bool ene_ok   = std::abs(result.get_energy()) <= std::abs(initial_mps.get_energy());
                        bool var_ok   = result.get_variance() <= initial_mps.get_variance();
                        bool rnorm_ok = result.get_rnorm() <= 1e-10;
                        // bool rnorm_nice = result.get_rnorm() <= 1e-12;
                        // if(rnorm_nice) break;
                        if(rnorm_ok and var_ok) break;

                        // The result is not good enough. But we could use it as initial guess
                        meta.optCost = OptCost::VARIANCE;
                        result       = internal::optimize_variance_eigs(tensors, result, status, meta);
                    }
                    break;
                }
                case OptAlgo::SUBSPACE: result = internal::optimize_variance_subspace(tensors, initial_mps, status, meta); break;
                case OptAlgo::SHIFTINV: throw except::runtime_error("OptCost::VARIANCE is not compatible with OptAlgo::SHIFTINV");
                case OptAlgo::MPSEIGS: throw except::runtime_error("OptAlgo::MPSEIGS has not been implemented yet");
                default: throw except::logic_error("Unrecognized OptAlgo::{}", static_cast<std::underlying_type_t<OptAlgo>>(meta.optAlgo));
            }
            break;
        }
    }

    tools::finite::opt::reports::print_eigs_report();

    // Finish up
    result.set_optsolver(meta.optSolver);
    result.set_optcost(meta.optCost);
    result.set_optalgo(meta.optAlgo);
    result.set_optexit(meta.optExit);
    result.validate_result();
    return result;
}

tools::finite::opt::opt_mps tools::finite::opt::find_ground_state(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                                  OptMeta &meta) {
    auto t_opt  = tid::tic_scope("opt");
    auto result = internal::optimize_energy_eigs(tensors, initial_mps, status, meta);
    tools::finite::opt::reports::print_eigs_report();
    // Finish up
    result.set_optsolver(meta.optSolver);
    result.set_optcost(meta.optCost);
    result.set_optalgo(meta.optAlgo);
    result.set_optexit(meta.optExit);
    result.validate_result();
    return result;
}

double tools::finite::opt::internal::windowed_func_abs(double x, double window) {
    if(std::abs(x) >= window)
        return std::abs(x) - window;
    else
        return 0;
}
double tools::finite::opt::internal::windowed_grad_abs(double x, double window) {
    if(std::abs(x) >= window)
        return num::sign(x);
    else
        return 0.0;
}

double tools::finite::opt::internal::windowed_func_pow(double x, double window) {
    if(std::abs(x) >= window)
        return x * x - window * window;
    else
        return 0.0;
}
double tools::finite::opt::internal::windowed_grad_pow(double x, double window) {
    if(std::abs(x) >= window)
        return 2.0 * x;
    else
        return 0.0;
}

std::pair<double, double> tools::finite::opt::internal::windowed_func_grad(double x, double window) {
    double func = 0;
    double grad = 0;
    if(std::abs(x) >= window) {
        if(x > 0) {
            func = (x - window) * (x - window);
            grad = 2 * (x - window);
        } else {
            func = (x + window) * (x + window);
            grad = 2 * (x + window);
        }
    }

    return std::make_pair(func, grad);
}

long tools::finite::opt::internal::get_ops(long d, long chiL, long chiR, long m) {
    if(chiR > chiL) return get_ops(d, chiR, chiL, m);
    if(d > m) {
        // d first
        long step1 = chiL * chiL * chiR * m * m * d;
        long step2 = chiL * chiR * d * m * m * m * (d * m + 1);
        long step3 = chiL * chiR * d * m * (chiR * m * m * m + m * m + 1);
        return step1 + step2 + step3;
    } else {
        // m first
        long step1 = chiL * chiL * chiR * m * d;
        long step2 = chiL * chiR * d * m * m * (d * m + 1);
        long step3 = chiL * chiR * chiR * d * m * (m + 1);
        return step1 + step2 + step3;
    }
}

namespace tools::finite::opt::internal {
    bool comparator::energy(const opt_mps &lhs, const opt_mps &rhs) {
        // The eigenvalue solver on H gives results sorted in energy
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_energy() < rhs.get_energy();
    }

    bool comparator::energy_distance(const opt_mps &lhs, const opt_mps &rhs, double target) {
        if(std::isnan(target)) throw except::logic_error("Energy target for comparison is NAN");
        auto diff_energy_lhs = std::abs(lhs.get_energy() - target);
        auto diff_energy_rhs = std::abs(rhs.get_energy() - target);
        return diff_energy_lhs < diff_energy_rhs;
    }

    bool comparator::variance(const opt_mps &lhs, const opt_mps &rhs) {
        // The eigenvalue solver on (H-E)² gives results sorted in variance
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_variance() < rhs.get_variance();
    }

    bool comparator::gradient(const opt_mps &lhs, const opt_mps &rhs) {
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_grad_max() < rhs.get_grad_max();
    }
    bool comparator::eigval(const opt_mps &lhs, const opt_mps &rhs) {
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_eigs_eigval() < rhs.get_eigs_eigval();
    }

    bool comparator::overlap(const opt_mps &lhs, const opt_mps &rhs) { return lhs.get_overlap() > rhs.get_overlap(); }

    bool comparator::eigval_and_overlap(const opt_mps &lhs, const opt_mps &rhs) {
        double ratio = std::max(lhs.get_eigs_eigval(), rhs.get_eigs_eigval()) / std::min(lhs.get_eigs_eigval(), rhs.get_eigs_eigval());
        if(ratio < 10 and lhs.get_overlap() >= std::sqrt(0.5)) return comparator::overlap(lhs, rhs);
        return comparator::eigval(lhs, rhs);
    }

    Comparator::Comparator(const OptMeta &meta_, double target_energy_) : meta(&meta_), target_energy(target_energy_) {}

    bool Comparator::operator()(const opt_mps &lhs, const opt_mps &rhs) {
        if(not meta) throw except::logic_error("No opt_meta given to comparator");
        // The variance case first
        if(meta->optCost == OptCost::VARIANCE) {
            return comparator::eigval_and_overlap(lhs, rhs);
        } else {
            auto diff = std::abs(lhs.get_eshift_eigval() - rhs.get_eshift_eigval());
            if(diff < settings::precision::eigs_tol_min) return lhs.get_overlap() > rhs.get_overlap();
            switch(meta->optRitz) {
                case OptRitz::NONE: throw std::logic_error("optimize_energy_eig_executor: Invalid: OptRitz::NONE");
                case OptRitz::SR: return comparator::energy(lhs, rhs);
                case OptRitz::LR: return comparator::energy(rhs, lhs);
                case OptRitz::SM: return comparator::energy_distance(lhs, rhs, 0.0);
                case OptRitz::IS: return comparator::energy_distance(lhs, rhs, target_energy);
                case OptRitz::TE: return comparator::energy_distance(lhs, rhs, target_energy);
            }
        }
        return comparator::eigval(lhs, rhs);
    }

    EigIdxComparator::EigIdxComparator(OptRitz ritz_, double shift_, double *data_, long size_) : ritz(ritz_), shift(shift_), eigvals(data_, size_) {}
    bool              EigIdxComparator::operator()(long lidx, long ridx) {
        auto lhs = eigvals[lidx];
        auto rhs = eigvals[ridx];
        switch(ritz) {
            case OptRitz::NONE: throw std::logic_error("EigvalComparator: Invalid OptRitz::NONE");
            case OptRitz::SR: return lhs < rhs;
            case OptRitz::LR: return rhs < lhs;
            case OptRitz::SM: return std::abs(lhs) < std::abs(rhs);
            case OptRitz::IS:
            case OptRitz::TE: {
                auto diff_lhs = std::abs(lhs - shift);
                auto diff_rhs = std::abs(rhs - shift);
                return diff_lhs < diff_rhs;
            }
            default: throw std::logic_error("EigvalComparator: Invalid OptRitz");
        }
    }
}
