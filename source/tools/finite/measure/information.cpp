
#include "../measure.h"
#include "config/settings.h"
#include "debug/info.h"
#include "general/iter.h"
#include "io/fmt_custom.h"
#include "math/eig.h"
#include "math/float.h"
#include "math/linalg/matrix.h"
#include "math/linalg/tensor.h"
#include "math/num.h"
#include "math/stat.h"
#include "math/tenx.h"
#include "qm/mpo.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/common/split.h"
#include "tools/finite/mps.h"
#include <Eigen/src/Core/GlobalFunctions.h>

template<typename ContainerType>
double count_finite(const ContainerType &arr) {
    double count = 0;
    for(const auto &a : arr) {
        if(std::isnan(a)) continue;
        if(std::isinf(a)) continue;
        count += 1;
    }
    return count;
}

template<typename ContainerType>
double nansum(const ContainerType &arr) {
    if(arr.size() == 0) return std::numeric_limits<double>::quiet_NaN();
    double sum   = 0;
    double count = 0;
    for(const auto &a : arr) {
        if(std::isnan(a)) continue;
        sum += a;
        count += 1;
    }
    if(count == 0) return std::numeric_limits<double>::quiet_NaN();
    return sum;
}

template<typename ContainerType>
double nanmean(const ContainerType &arr) {
    if(arr.size() == 0) return std::numeric_limits<double>::quiet_NaN();
    double sum   = 0;
    double count = 0;
    for(const auto &a : arr) {
        if(std::isnan(a)) continue;
        sum += a;
        count += 1;
    }
    if(count == 0) return std::numeric_limits<double>::quiet_NaN();
    return sum / count;
}

double tools::finite::measure::entanglement_entropy_log2(const StateFinite &state, size_t nsites /* sites to the left of the partition */) {
    auto t_ent = tid::tic_scope("neumann_entropy", tid::level::highest);
    if(nsites == 0ul) return 0.0;
    auto pos     = state.get_position<long>();
    auto pos_tgt = 0ul;
    if(pos + 1ul < nsites) {
        // Get the L of the B-site at pos == nsites-1
        pos_tgt = nsites - 1ul;
    } else if(pos + 1ul == nsites) {
        // Get the LC of the AC site at pos == nsites-1
        pos_tgt = nsites - 1ul;
    } else {
        // Get the L of the A (or AC) -site at pos == nsites
        pos_tgt = nsites;
    }
    const auto &mps = state.get_mps_site(pos_tgt);
    if(mps.isCenter()) {
        auto &LC = mps.get_LC();
        auto  SE = Eigen::Tensor<cplx, 0>(-LC.square().contract(LC.square().log2().eval(), tenx::idx({0}, {0})));
        return std::abs(SE(0));
    } else {
        auto &L  = mps.get_L();
        auto  SE = Eigen::Tensor<cplx, 0>(-L.square().contract(L.square().log2().eval(), tenx::idx({0}, {0})));
        return std::abs(SE(0));
    }
}

std::vector<double> tools::finite::measure::entanglement_entropies_log2(const StateFinite &state) {
    auto                t_ent = tid::tic_scope("neumann_entropy", tid::level::highest);
    std::vector<double> entanglement_entropies;
    entanglement_entropies.reserve(state.get_length() + 1);
    if(not state.has_center_point()) entanglement_entropies.emplace_back(0);
    for(const auto &mps : state.mps_sites) {
        auto                  &L  = mps->get_L();
        Eigen::Tensor<cplx, 0> SE = -L.square().contract(L.square().log2().eval(), tenx::idx({0}, {0}));
        entanglement_entropies.emplace_back(std::abs(SE(0)));
        if(mps->isCenter()) {
            auto &LC = mps->get_LC();
            SE       = -LC.square().contract(LC.square().log2().eval(), tenx::idx({0}, {0}));
            entanglement_entropies.emplace_back(std::abs(SE(0)));
        }
    }
    if(entanglement_entropies.size() != state.get_length() + 1) throw except::logic_error("entanglement_entropies.size() should be length+1");
    if(entanglement_entropies.front() != 0.0) throw except::logic_error("First entropy should be 0. Got: {:.16f}", entanglement_entropies.front());
    if(entanglement_entropies.back() != 0.0) throw except::logic_error("Last entropy should be 0. Got: {:.16f}", entanglement_entropies.back());
    return entanglement_entropies;
}

double tools::finite::measure::subsystem_entanglement_entropy_log2(const StateFinite &state, const std::vector<size_t> &sites) {
    if(sites.empty()) return 0;
    // If sites is contiguous, then we can return a bipartite entropy
    if(sites == num::range<size_t>(0ul, sites.size())) { return tools::finite::measure::entanglement_entropy_log2(state, sites.size()); }
    bool is_real = state.is_real();
    auto solver  = eig::solver();
    auto evs     = Eigen::ArrayXd();
    if(is_real) {
        auto rho = state.get_reduced_density_matrix<real>(sites);
        tools::log->trace("eig rho_real: {} ...", rho.dimensions());
        if(rho.dimension(0) <= 8192) {
            solver.eig<eig::Form::SYMM>(rho.data(), rho.dimension(0), eig::Vecs::OFF);
        } else {
            solver.eig<eig::Form::SYMM>(rho.data(), rho.dimension(0), 'V', -1, -1, 1e-14, 1, eig::Vecs::OFF);
        }
        evs = eig::view::get_eigvals<real>(solver.result); // Eigenvalues of rho
        tools::log->trace("eig rho_real: {} ... time {:.2e} s (prep {:.2e} s)", rho.dimensions(), solver.result.meta.time_total, solver.result.meta.time_prep);
    } else {
        auto rho = state.get_reduced_density_matrix<cplx>(sites);
        tools::log->trace("eig rho_cplx: {} ...", rho.dimensions());

        solver.eig<eig::Form::SYMM>(rho.data(), rho.dimension(0), eig::Vecs::OFF);
        evs = eig::view::get_eigvals<real>(solver.result).real(); // Eigenvalues of rho
        tools::log->trace("eig rho_cplx: {} ... time {:.2e} s (prep {:.2e} s)", rho.dimensions(), solver.result.meta.time_total, solver.result.meta.time_prep);
    }
    double entanglement_entropy_log2 = 0;
    for(const auto &e : evs) {
        if(e > 0) entanglement_entropy_log2 += -e * std::log2(e); // We use log base 2 for information
    }
    return entanglement_entropy_log2;
}

Eigen::ArrayXXd tools::finite::measure::subsystem_entanglement_entropies_log2(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    if(state.get_length<long>() >= 10) {
        return subsystem_entanglement_entropies_swap_log2(state, svd_cfg, uc); // trnc lim 1e-6 is 12 digits correct on SE
    }
    if(uc == UseCache::TRUE and state.measurements.subsystem_entanglement_entropies.has_value())
        return state.measurements.subsystem_entanglement_entropies.value();
    auto            t_see  = tid::tic_scope("subsystem_entanglement_entropies", tid::level::normal);
    auto            len    = state.get_length<long>();
    auto            bee    = measure::entanglement_entropies_log2(state); // bipartite entanglement entropies
    auto            solver = eig::solver();
    Eigen::ArrayXXd see    = Eigen::ArrayXXd::Zero(len, len); // susbsystem entanglement entropy
    for(long off = 0; off < len; ++off) {
        for(long ext = 1; ext <= len; ++ext) {
            if(off + ext > len) continue;
            // Check if the segment includes the edge, in which case
            // we can simply take the bipartite entanglement entropy
            bool is_left_edge  = off == 0;
            bool is_right_edge = off + ext == len;
            // tools::log->trace("Calculating subsystem entanglement entropy: extent (l): {} | offset (n) {}", ext, off);
            if(is_left_edge) {
                auto idx          = safe_cast<size_t>(off + ext);
                see(ext - 1, off) = bee[idx];
            } else if(is_right_edge) {
                auto idx          = safe_cast<size_t>(off);
                see(ext - 1, off) = bee[idx];
            } else {
                auto sites        = num::range<size_t>(off, off + ext);
                auto evs          = Eigen::ArrayXd();
                see(ext - 1, off) = subsystem_entanglement_entropy_log2(state, sites);
            }
            tools::log->debug("RSS {:.3f} MB | Peak {:.3f} MB", debug::mem_rss_in_mb(), debug::mem_hwm_in_mb());
        }
        // if(uc == UseCache::TRUE) state.clear_cache();
    }
    // tools::log->info("subsystem entanglement entropies {:.3e} s: \n{}\n", t_see->get_last_interval(), linalg::matrix::to_string(see, 16));

    if(uc == UseCache::TRUE) state.measurements.subsystem_entanglement_entropies = see;
    return see;
}

Eigen::ArrayXXd tools::finite::measure::subsystem_entanglement_entropies_swap_log2(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    if(uc == UseCache::TRUE and state.measurements.subsystem_entanglement_entropies.has_value())
        return state.measurements.subsystem_entanglement_entropies.value();
    svd_cfg                   = svd_cfg.value_or(svd::config());
    svd_cfg->svd_lib          = svd_cfg->svd_lib.value_or(svd::lib::lapacke);
    svd_cfg->svd_rtn          = svd_cfg->svd_rtn.value_or(svd::rtn::gesdd);
    svd_cfg->rank_max         = svd_cfg->rank_max.value_or(settings::storage::dataset::subsystem_entanglement_entropies::bond_lim);
    svd_cfg->truncation_limit = svd_cfg->truncation_limit.value_or(settings::storage::dataset::subsystem_entanglement_entropies::trnc_lim);
    auto            t_see     = tid::tic_scope("subsystem_entanglement_entropies", tid::level::normal);
    auto            len       = state.get_length<long>();
    auto            solver    = eig::solver();
    Eigen::ArrayXXd see       = Eigen::ArrayXXd::Zero(len, len); // susbsystem entanglement entropy

    auto state_swap = state;
    auto length     = state_swap.get_length<size_t>();
    auto sites      = num::range<size_t>(0, length);
    // Truncate the state down to the bond and truncation limits asked for first
    tools::log->info("subsystem_entanglement_entropies_swap_log2: Normalizing state | svd: {}", svd_cfg->to_string());
    finite::mps::normalize_state(state_swap, svd_cfg, NormPolicy::ALWAYS);
    finite::mps::move_center_point_to_pos_dir(state_swap, 0, 1, svd_cfg);
    tools::log->info("calculating subsystem entanglement entropies with the swap method | svd: {}", svd_cfg->to_string());
    for(long off = 0; off < len; ++off) {
        // tools::log->debug("subsystem_entanglement_entropies_swap_log2: calculating offset {} ...", off);
        if(off > 0) {
            // Move the left-most site to one site beyond the last site
            // Example with 6 sites
            //      012345 --> 123450
            //      123450 --> 234510
            for(size_t i = 0; i < length - safe_cast<size_t>(off); ++i) {
                // auto [bond_old_L, bond_old_R] = measure::bond_dimensions(state_swap, i);
                // auto trnc_old                 = measure::truncation_errors(state_swap);
                // auto sval_old                 = state_swap.get_bond(i, i + 1);
                mps::swap_sites(state_swap, i, i + 1, sites, GateMove::OFF, svd_cfg);
                // auto [bond_new_L, bond_new_R] = measure::bond_dimensions(state_swap, i);
                // auto trnc_new                 = measure::truncation_errors(state_swap);
                // auto sval_new                 = state_swap.get_bond(i, i + 1);

                // tools::log->info("swap off {:2} | pos {} |  {} <-> {} | sites {::2} | bonds {::3}", off, state_swap.get_position<long>(), i, i + 1, sites,
                //                  measure::bond_dimensions(state_swap));
                // if(bond_new_R != sval_new.size())
                // throw except::logic_error("bond dimension mismatch: bond_new_R {} != sval_new.size() {}", bond_new_R, sval_new.size());
                // for(long sidx = 0; sidx < std::max(sval_old.size(), sval_new.size()); ++sidx) {
                // auto sold = sidx < sval_old.size() ? sval_old.coeff(sidx).real() : 0.0;
                // auto snew = sidx < sval_new.size() ? sval_new.coeff(sidx).real() : 0.0;
                // tools::log->info("idx {:>4}: {:>10.5e} --> {:>10.5e}", sidx, sold, snew);
                // }
            }
        }
        tools::log->info("swap off {:2} | pos {} | sites {::2} | bonds {::3}", off, state_swap.get_position<long>(), sites,
                         measure::bond_dimensions(state_swap));
        auto bee = measure::entanglement_entropies_log2(state_swap); // bipartite entanglement entropies

        for(long ext = 1; ext <= len; ++ext) {
            if(off + ext > len) continue;
            auto idx          = safe_cast<size_t>(ext);
            see(ext - 1, off) = bee[idx];
        }
    }
    // tools::log->info("subsystem entanglement entropies swap {:.3e} s: \n{}\n", t_see->get_last_interval(), linalg::matrix::to_string(see, 16));
    // tools::log->debug("RSS {:.3f} MB | Peak {:.3f} MB", debug::mem_rss_in_mb(), debug::mem_hwm_in_mb());
    if(uc == UseCache::TRUE) state.measurements.subsystem_entanglement_entropies = see;
    return see;
}

/*! Calculates the subsystem entanglement entropies,
 * starting with the cheapest entries and progressively harder until a threshold number of bits have been found.
 * Use max_bit_error to set the threshold:
 *  - use a positive value to indicate relative error
 *  - use a negative value to indicate absolute error (i.e. how many missing bits you can accept)
 */
Eigen::ArrayXXd tools::finite::measure::subsystem_entanglement_entropies_fast_log2(const StateFinite &state, double max_bit_error,
                                                                                   std::optional<svd::config> svd_cfg, UseCache uc) {
    if(uc == UseCache::TRUE and state.measurements.subsystem_entanglement_entropies.has_value())
        return state.measurements.subsystem_entanglement_entropies.value();
    svd_cfg                   = svd_cfg.value_or(svd::config());
    svd_cfg->svd_lib          = svd_cfg->svd_lib.value_or(svd::lib::lapacke);
    svd_cfg->svd_rtn          = svd_cfg->svd_rtn.value_or(svd::rtn::gesdd);
    svd_cfg->rank_max         = svd_cfg->rank_max.value_or(settings::storage::dataset::subsystem_entanglement_entropies::bond_lim);
    svd_cfg->truncation_limit = svd_cfg->truncation_limit.value_or(settings::storage::dataset::subsystem_entanglement_entropies::trnc_lim);
    auto t_see                = tid::tic_scope("subsystem_entanglement_entropies", tid::level::normal);
    auto state_temp           = state;
    auto length               = state_temp.get_length<size_t>();
    auto sites                = num::range<size_t>(0, length);
    // Truncate the state down to the bond and truncation limits asked for first
    tools::log->info("subsystem_entanglement_entropies_swap_log2: Normalizing state | svd: {}", svd_cfg->to_string());
    finite::mps::normalize_state(state_temp, svd_cfg, NormPolicy::ALWAYS);
    finite::mps::move_center_point_to_pos_dir(state_temp, 0, 1, svd_cfg);
    tools::log->info("calculating subsystem entanglement entropies with the swap method | svd: {}", svd_cfg->to_string());
    constexpr auto  nan = std::numeric_limits<double>::quiet_NaN();
    auto            len = state_temp.get_length<long>();
    auto            bee = measure::entanglement_entropies_log2(state_temp); // bipartite entanglement entropies
    Eigen::ArrayXXd see = Eigen::ArrayXXd::Zero(len, len) * nan;            // susbsystem entanglement entropies

    // We begin with the low-hanging fruit (1): Bipartite entanglement entropies
    for(long ext = 1; ext <= len; ++ext) { // Assume off == 0
        auto idx                = safe_cast<size_t>(ext);
        see(ext - 1, 0)         = bee[idx];
        see(len - ext - 1, ext) = bee[idx];
    }
    // Next we take small subsystems
    long max_ext = 11; // max rho size 2048x2048
    for(long off = 1; off < len; ++off) {
        for(long ext = 1; ext <= max_ext; ++ext) {
            if(off + ext >= len) continue; // We can have equality because we get the diagonal in see for free
            auto subsystem    = num::range<size_t>(off, off + ext);
            see(ext - 1, off) = subsystem_entanglement_entropy_log2(state_temp, subsystem);
        }
    }

    // Finally, we use the swap algorithm for the missing entries until we find most bits
    for(long off = 1; off < len - max_ext - 1; ++off) { // -1 because the diagonal in see is equal to its left-most column.
        // if(off + max_ext > len) continue;

        // Update the information lattice
        auto info       = information_lattice(see);
        auto bits_total = state_temp.get_length<double>();
        auto bits_found = info.isNaN().select(0, info).sum();
        auto bits_error = 1 - bits_found / bits_total;

        tools::log->info("see: \n{}\n", linalg::matrix::to_string(see, 10));
        tools::log->info("info: \n{}\n", linalg::matrix::to_string(info, 10));
        tools::log->info("bits found {:.16f}/{:.16f} | error {:.16f} | com {:.16f}", bits_found, bits_total, bits_error, information_center_of_mass(info));
        if(max_bit_error >= 0 and bits_error <= max_bit_error) break;
        if(max_bit_error < 0 and std::abs(bits_total - bits_found) <= std::abs(max_bit_error)) break;

        // tools::log->debug("subsystem_entanglement_entropies_swap_log2: calculating offset {} ...", off);
        // Move the left-most site to one site beyond the last site
        // Example with 6 sites
        //      012345 --> 123450
        //      123450 --> 234510
        for(size_t i = 0; i < length - safe_cast<size_t>(off); ++i) {
            // auto [bond_old_L, bond_old_R] = measure::bond_dimensions(state_swap, i);
            // auto trnc_old                 = measure::truncation_errors(state_swap);
            // auto sval_old                 = state_swap.get_bond(i, i + 1);
            mps::swap_sites(state_temp, i, i + 1, sites, GateMove::OFF, svd_cfg);
        }
        tools::log->info("swap off {:2} | pos {} | sites {::2} | bonds {::3}", off, state_temp.get_position<long>(), sites,
                         measure::bond_dimensions(state_temp));
        auto bee_swap = measure::entanglement_entropies_log2(state_temp); // bipartite entanglement entropies for the swapped state
        for(long ext = max_ext + 1; ext <= len; ++ext) {
            if(off + ext > len) continue;
            auto idx          = safe_cast<size_t>(ext);
            see(ext - 1, off) = bee_swap[idx];
        }
    }

    if(uc == UseCache::TRUE) state.measurements.subsystem_entanglement_entropies = see;
    return see;
}

Eigen::ArrayXXd tools::finite::measure::information_lattice(const Eigen::ArrayXXd &SEE) {
    auto            Ldim        = SEE.rows();
    Eigen::ArrayXXd infolattice = Eigen::ArrayXXd::Zero(Ldim, Ldim);
    for(long nidx = 0; nidx < Ldim; ++nidx) /* The offset */ {
        for(long lidx = 0; lidx < Ldim; ++lidx) /* The extent, or "scale" */ {
            if(lidx + nidx >= Ldim) continue;
            if(lidx == 0) {
                infolattice(lidx, nidx) = 1 - SEE(lidx, nidx);
            } else if(lidx == 1) {
                infolattice(lidx, nidx) = -SEE(lidx, nidx) + SEE(lidx - 1, nidx) + SEE(lidx - 1, nidx + 1);
            } else {
                infolattice(lidx, nidx) = -SEE(lidx, nidx) + SEE(lidx - 1, nidx) + SEE(lidx - 1, nidx + 1) - SEE(lidx - 2, nidx + 1);
            }
        }
    }
    return infolattice;
}

Eigen::ArrayXXd tools::finite::measure::information_lattice(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    if(uc == UseCache::TRUE and state.measurements.information_lattice.has_value()) return state.measurements.information_lattice.value();
    auto SEE         = subsystem_entanglement_entropies_log2(state, svd_cfg, uc);
    auto infolattice = information_lattice(SEE);
    if(uc == UseCache::TRUE) state.measurements.information_lattice = infolattice;
    return infolattice;
}

Eigen::ArrayXd tools::finite::measure::information_per_scale(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    if(uc == UseCache::TRUE and state.measurements.information_per_scale.has_value()) return state.measurements.information_per_scale.value();
    Eigen::ArrayXd info_per_scale = information_lattice(state, svd_cfg, uc).rowwise().sum();
    if(uc == UseCache::TRUE) state.measurements.information_per_scale = info_per_scale;
    // tools::log->info("info per scale: {::.3e}", state.measurements.information_per_scale.value());
    return info_per_scale;
}

/*! For masses m_i at coordinates x_i, the center of mass is defined as
 *      sum_i m_i * x_i
 *      ---------------
 *      sum_i m_i
 */
double tools::finite::measure::information_center_of_mass(const Eigen::ArrayXXd &information_lattice) {
    Eigen::ArrayXd il = information_lattice.isNaN().select(0, information_lattice).rowwise().sum();
    auto           l  = Eigen::ArrayXd::LinSpaced(il.size(), 0.0, static_cast<double>(il.size() - 1l));
    auto           ic = il.cwiseProduct(l).sum() / il.sum();
    return ic;
}
/*! For masses m_i at coordinates x_i, the center of mass is defined as
 *      sum_i m_i * x_i
 *      ---------------
 *      sum_i m_i
 */
double tools::finite::measure::information_center_of_mass(const Eigen::ArrayXd &information_per_scale) {
    const auto &il = information_per_scale;
    auto        l  = Eigen::ArrayXd::LinSpaced(il.size(), 0.0, static_cast<double>(il.size() - 1l));
    auto        ic = il.cwiseProduct(l).sum() / il.sum();
    return ic;
}

/*! For masses m_i at coordinates x_i, the center of mass is defined as
 *      sum_i m_i * x_i
 *      ---------------
 *      sum_i m_i
 */
double tools::finite::measure::information_center_of_mass(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    return information_center_of_mass(information_per_scale(state, svd_cfg, uc));
}

/*! Calculates xi from the expectation value of the geometric distribution, such that p = 1-exp(-1/xi) */
double tools::finite::measure::information_xi_from_geometric_dist(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    Eigen::ArrayXd il = information_per_scale(state, svd_cfg, uc);
    auto           L  = il.size();
    if(il.size() <= 2) return std::numeric_limits<double>::quiet_NaN();
    Eigen::ArrayXd l    = Eigen::ArrayXd::LinSpaced(L, 0.0, static_cast<double>(L - 1l));
    Eigen::ArrayXd mask = Eigen::ArrayXd::Ones(L); // Mask to select the interesting interval

    bool has_extended_info = il.bottomRows(L / 2).sum() >= 0.9; // Both topological and thermal states have >= 1 bits beyond L/2.
    bool decay_starts_at_0 = il[0] > il[1];
    bool decay_starts_at_1 = il[1] >= il[0] and il[1] > il[2];

    if(has_extended_info) mask.bottomRows(L / 2) = 0; // If there is a topological bit at l ~ L? we should only consider l in [0,L/2]
    if(decay_starts_at_1) mask[0] = 0;                // If the decay starts at l==1, take the geometric distribution starting at l==1

    auto ic = (il * mask).cwiseProduct(l).sum() / (il * mask).sum(); // information center of mass (or <l>)
    auto p  = 1.0;                                                   // Probability p for the geometric distribution

    if(decay_starts_at_0) {
        p = 1 / (1 + ic); // The geometric decay starts at l == 0
    } else if(decay_starts_at_1) {
        p = 1 / ic; // The geometric decay starts at l == 1
    }
    auto xi = -1.0 / std::log(1 - p);
    // auto N  = decay_starts_at_1 ? il[1] : il[0]; // The initial value of the exponential decay
    // auto lstart = decay_starts_at_1 ? 1.0 : 0.0;
    // auto y_geom = Eigen::ArrayXd((il * mask).sum() * p * pow(1 - p, l - lstart));
    // auto y_expd = Eigen::ArrayXd(N * mask * exp(-(l - lstart) / xi));
    // tools::log->info("info_xi_from_geom: ic       : {:.16f}", ic);
    // tools::log->info("info_xi_from_geom: p        : {:.16f}", p);
    // tools::log->info("info_xi_from_geom: xi       : {:.16f}", xi);
    // tools::log->info("info_xi_from_geom: y_info   : {::.16f}", il);
    // tools::log->info("info_xi_from_geom: y_geom   : {::.16f}", y_geom);
    // tools::log->info("info_xi_from_geom: y_expfit : {::.16f}", y_expd);
    return xi;
}

/*! Calculates xi from the slope of the decay in the exponential distribution in log-scale */
double tools::finite::measure::information_xi_from_avg_log_slope(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::ArrayXd il  = information_per_scale(state, svd_cfg, uc);
    if(il.size() <= 3) return nan;
    auto           L    = il.size();
    Eigen::ArrayXd l    = Eigen::ArrayXd::LinSpaced(L, 0.0, static_cast<double>(L - 1l));
    Eigen::ArrayXd mask = Eigen::ArrayXd::Ones(L); // Mask to select the interesting interval

    bool has_extended_info = il.bottomRows(L / 2).sum() >= 0.9; // Both topological and thermal states have >= 1 bits beyond L/2.
    bool decay_starts_at_1 = il[1] >= il[0] and il[1] > il[2];

    if(has_extended_info) mask.bottomRows(L / 2 - 1) = nan; // If there is a topological bit at l ~ L? we should only consider l in [0,L/2]
    if(decay_starts_at_1) mask[0] = nan;                    // If the decay starts at l==1, take the geometric distribution starting at l==1
    mask = (il < 1e-15).select(nan, mask);                  // Filter tiny noisy values

    auto il_log      = Eigen::ArrayXd(Eigen::log(il * mask));
    auto il_log_diff = Eigen::ArrayXd(il_log.bottomRows(L - 1) - il_log.topRows(L - 1));
    auto xi          = -1.0 / nanmean(il_log_diff);

    // auto lstart = decay_starts_at_1 ? 1.0 : 0.0;
    // auto istart = decay_starts_at_1 ? 1 : 0;
    // auto y_expd = Eigen::ArrayXd(il[istart] * mask * exp(-(l - lstart) / xi));
    // tools::log->info("info_xi_from_expd: xi            : {:.16f}", xi);
    // tools::log->info("info_xi_from_expd: y_info        : {::.16f}", il);
    // tools::log->info("info_xi_from_expd: y_expfit      : {::.16f}", y_expd);
    // tools::log->info("info_xi_from_expd: il_log        : {::.6f}", il_log);
    // tools::log->info("info_xi_from_expd: il_log_diff   : {::.6f}", il_log_diff);
    // tools::log->info("info_xi_from_expd: il_log_diff⁻¹ : {::.6f}", il_log_diff.cwiseInverse());
    return xi;
}

double tools::finite::measure::information_xi_from_exp_fit(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::ArrayXd il  = information_per_scale(state, svd_cfg, uc);
    if(il.size() <= 3) return nan;
    auto           L    = il.size();
    Eigen::ArrayXd l    = Eigen::ArrayXd::LinSpaced(L, 0.0, static_cast<double>(L - 1l));
    Eigen::ArrayXd mask = Eigen::ArrayXd::Ones(L); // Mask to select the interesting interval

    bool has_extended_info = il.bottomRows(L / 2).sum() >= 0.9; // Both topological and thermal states have >= 1 bits beyond L/2.
    bool decay_starts_at_1 = il[1] >= il[0] and il[1] > il[2];

    if(has_extended_info) mask.bottomRows(L / 2 - 1) = 0; // If there is a topological bit at l ~ L? we should only consider l in [0,L/2]
    if(decay_starts_at_1) mask[0] = 0;                    // If the decay starts at l==1, take the geometric distribution starting at l==1
    mask = (il < 1e-15).select(0, mask);                  // Filter tiny noisy values

    // Collect the well-defined points in the domain
    auto Y = std::vector<double>();
    auto X = std::vector<double>();
    for(const auto &[i, m] : iter::enumerate(mask))
        if(m == 1) {
            Y.emplace_back(std::log(il[i]));
            X.emplace_back(l[i]);
        }
    tools::log->info("X : {}", X);
    tools::log->info("Y : {::.16f}", Y);
    auto linfit = stat::linearFit(X, Y);
    return -1.0 / linfit.slope;
}

double tools::finite::measure::information_scale_by_frac(const StateFinite &state, double fraction, std::optional<svd::config> svd_cfg, UseCache uc) {
    // fraction = 0.5 implies the scale at which we find half of the total information
    assert(fraction >= 0.0);
    assert(fraction <= 1.0);
    auto info_per_scale = information_per_scale(state, svd_cfg, uc);
    auto L              = info_per_scale.size();
    auto idx_interp     = num::LinSpaced(1000 * L, 0, L);
    auto bit_interp     = num::interp1_prev(num::range(0, L), info_per_scale, idx_interp);
    auto sum_interp     = num::cumtrapz(idx_interp, bit_interp);

    double bits_to_count  = fraction * static_cast<double>(L);
    size_t idx_infocenter = 0;
    for(long idx = 0; idx < sum_interp.size(); ++idx) {
        if(sum_interp[idx] > bits_to_count) break;
        idx_infocenter = idx;
    }
    return idx_interp[idx_infocenter];
}

double tools::finite::measure::information_typ_scale(const StateFinite &state, std::optional<svd::config> svd_cfg, UseCache uc) {
    if(uc == UseCache::TRUE and state.measurements.information_typ_scale.has_value()) return state.measurements.information_typ_scale.value();
    auto info_typ_scale = information_scale_by_frac(state, 0.5, svd_cfg, uc) * std::log(2);
    if(uc == UseCache::TRUE) state.measurements.information_typ_scale = info_typ_scale;
    return info_typ_scale;
}
