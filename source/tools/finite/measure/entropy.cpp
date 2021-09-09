
#include "../measure.h"
#include <bitset>
#include <config/debug.h>
#include <general/iter.h>
#include <math/tenx.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/mps.h>
#include <utility>

struct Amplitude {
    long                                  state_size;          // System size
    std::bitset<64>                       bits;                // Bits that select spins on each MPS site
    std::optional<long>                   site = std::nullopt; // MPS site (not Schmidt site!)
    Eigen::Tensor<StateFinite::Scalar, 1> ampl;                // Accumulates the MPS tensors
    Amplitude(long state_size, std::bitset<64> bits, std::optional<long> site, Eigen::Tensor<StateFinite::Scalar, 1> ampl)
        : state_size(state_size), bits(bits), site(site), ampl(std::move(ampl)) {}

    private:
    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_string(const std::bitset<64> &b, long num) const {
        if constexpr(on) {
            if(tools::log->level() > spdlog::level::trace) return {};
            assert(num > 0 and num < 64 and "num should be in range [0,64]");
            return fmt::format(FMT_STRING("{1:>{0}}"), state_size, b.to_string().substr(b.size() - num));
        } else
            return {};
    }
    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_rstring(const std::bitset<64> &b, long num) const {
        if constexpr(on) {
            if(tools::log->level() > spdlog::level::trace) return {};
            assert(num > 0 and num < 64 and "num should be in range [0,64]");
            auto bs = b.to_string().substr(b.size() - num);
            return fmt::format(FMT_STRING("{1:<{0}}"), state_size, std::string{bs.rbegin(), bs.rend()});
        } else
            return {};
    }

    public:
    [[nodiscard]] std::string to_string() const {
        if(site)
            return to_string(bits, state_size);
        else
            return {};
    }
    [[nodiscard]] std::string to_rstring() const {
        if(site)
            return to_rstring(bits, state_size);
        else
            return {};
    }
    [[nodiscard]] std::string to_string(bool reverse) const {
        if(site)
            if(reverse)
                return to_rstring(bits, state_size);
            else
                return to_string(bits, state_size);
        else
            return {};
    }

    void eval_from_A(const StateFinite &state, long tgt_pos, std::vector<Amplitude> &cache) {
        if(not site or site.value() < tgt_pos) {
            auto t_evalA = tid::tic_scope("evalA");

            // There are missing mps in the amplitude. Let's add them
            // First, check if it is already available in the amplitudes database
            // We can look it up by transforming the bits between 0 to bit_site to an integer
            // and use that as the index on amplitudes
            long state_pos = state.get_position<long>();
            if(tgt_pos > state_pos) throw std::logic_error(fmt::format("eval_from_A: expected mps_site ({}) <= state_pos ({})", tgt_pos, state_pos));

            if(tgt_pos > 0 and not cache.empty()) {
                // There is a chance to continue building an existing amplitude
                auto            t_check    = tid::tic_scope("cache_check");
                std::bitset<64> bits_index = 0;
                for(size_t i = 0; i < static_cast<size_t>(tgt_pos); i++) bits_index[i] = bits[i];
                if(cache.size() > bits_index.to_ulong()) {
                    const auto &c = cache.at(bits_index.to_ulong()); // A cache item
                    if(c.site and c.site.value() <= tgt_pos) {
                        // Cache hit! No need to compute the amplitude from scratch, just append the next site
                        auto t_hit = tid::tic_scope("cache_hit");
                        site       = c.site;
                        ampl       = c.ampl;
                        if constexpr(settings::debug_numen)
                            tools::log->trace(FMT_STRING("from A: cache hit: [site {:>2} | n {:>2} | bits {}] target [site {:>2} | n {:>2} | bits {}"),
                                              c.site.value(), bits_index.count(), to_rstring(bits_index, c.site.value() + 1), tgt_pos, bits.count(),
                                              to_rstring(bits, tgt_pos + 1));
                    }
                }
            }

            if(ampl.size() == 0) { // Initialize
                ampl.resize(1);
                ampl.setConstant(1.0);
            }
            Eigen::Tensor<StateFinite::Scalar, 1> temp;
            // Contract the missing mps up to, but not including, the last mps at mps_site
            for(const auto &mps : state.mps_sites) {
                long pos = mps->get_position<long>();
                if(site and site.value() >= pos) continue; // Fast-forward to the missing sites
                if(tgt_pos < pos) break;                   // Contract up to the mps at tgt_pos
                if(ampl.size() != mps->get_chiL())
                    throw std::runtime_error(fmt::format(FMT_STRING("eval() failed for site {}: "
                                                                    "mismatch in ampl({}) with size = {} and mps({}) with chiL = {} | bits {}"),
                                                         tgt_pos, site.value(), ampl.size(), pos, mps->get_chiL(), to_rstring()));

                long                size = mps->get_chiR();
                std::array<long, 3> off  = {bits[static_cast<size_t>(pos)], 0, 0}; // This selects which bit gets appended to ampl
                std::array<long, 3> ext  = {1, mps->get_chiL(), mps->get_chiR()};
                // ampl never has a trailing Lambda
                auto t_con = tid::tic_scope("contract");
                if constexpr(settings::debug_numen)
                    tools::log->trace(FMT_STRING("from A: cntrction: pos  {:>2} | tgt {:>2} | bits {}"), pos, tgt_pos, to_rstring(bits, pos));

                temp.resize(size);
                temp.device(tenx::omp::getDevice()) = ampl.contract(mps->get_M_bare().slice(off, ext), tenx::idx({0}, {1})).reshape(std::array<long, 1>{size});
                t_con.toc();
                ampl = temp;

                // Update the current site
                site = pos;
            }
            if(not site) throw std::logic_error(fmt::format("ampl has undefined site: should be {}", tgt_pos));
            if(site.value() != tgt_pos) throw std::logic_error(fmt::format("site ({}) != mps_site ({})", site.value(), tgt_pos));
        }
    }

    void eval_from_B(const StateFinite &state, long tgt_pos, std::vector<Amplitude> &cache) {
        // Start by calculating the mps sites that should be included in the amplitude
        // Remember that in eval_from_B we calculate the amplitude starting from the right-end of the chain
        if(not site or site.value() > tgt_pos) {
            auto t_evalB = tid::tic_scope("evalB");

            // There are missing mps in the amplitude. Let's add them
            // First, check if it is already available in the amplitudes database
            // We can look it up by transforming the bits between 0 to bit_site to an integer
            // and use that as the index on amplitudes
            long state_len = state.get_length<long>();
            long state_pos = state.get_position<long>();
            long mps_rsite = state_len - 1 - tgt_pos;
            if(tgt_pos <= state_pos) throw std::logic_error(fmt::format("eval_from_B: expected mps_site ({}) > state_pos ({})", tgt_pos, state_pos));

            if(tgt_pos < state_len - 1 and not cache.empty()) {
                // There is a chance to continue building an existing amplitude
                auto            t_check    = tid::tic_scope("cache_check");
                std::bitset<64> bits_index = 0;
                for(size_t i = 0; i < static_cast<size_t>(mps_rsite); i++) bits_index[i] = bits[i];
                if(cache.size() > bits_index.to_ulong()) {
                    const auto &c = cache.at(bits_index.to_ulong()); // A cache item
                    if(c.site and c.site.value() > tgt_pos) {
                        // Cache hit! No need to compute the amplitude from scratch
                        auto t_hit = tid::tic_scope("cache_hit");
                        site       = c.site;
                        ampl       = c.ampl;
                        if constexpr(settings::debug_numen)
                            tools::log->trace(FMT_STRING("from B: cache hit: [sites {:>2} | n {:>2} | bits {}] target [sites {:>2} | n {:>2} | bits {}"),
                                              c.site.value(), bits_index.count(), to_string(bits_index, c.site.value() + 1), tgt_pos, bits.count(),
                                              to_string(bits, tgt_pos + 1));
                    }
                }
            }

            if(ampl.size() == 0) { // Initialize
                ampl.resize(1);
                ampl.setConstant(1.0);
            }
            Eigen::Tensor<StateFinite::Scalar, 1> temp;
            // Contract the missing mps
            for(const auto &mps : iter::reverse(state.mps_sites)) {
                long pos = mps->get_position<long>();
                if(site and site.value() <= pos) continue; // Fast-forward to the missing sites
                if(tgt_pos > pos) break;                   // Contract up to the mps at mps_pos
                if(ampl.size() != mps->get_chiR())
                    throw std::runtime_error(fmt::format(FMT_STRING("eval() failed for site {}: "
                                                                    "mismatch in ampl({}) with size = {} and mps({}) with chiR = {} | bits {}"),
                                                         tgt_pos, site.value(), ampl.size(), pos, mps->get_chiL(), to_string()));
                long                mps_rpos = state_len - 1 - pos;
                long                size     = mps->get_chiL();
                std::array<long, 3> off      = {bits[static_cast<size_t>(mps_rpos)], 0, 0}; // This selects which bit gets prepended to ampl
                std::array<long, 3> ext      = {1, mps->get_chiL(), mps->get_chiR()};

                // ampl never has a trailing Lambda
                auto t_con = tid::tic_scope("contract");
                if constexpr(settings::debug_numen)
                    tools::log->trace(FMT_STRING("from B: cntrction: pos  {:>2} | tgt {:>2} | bits {}"), pos, tgt_pos, to_string(bits, state_len - pos));
                temp.resize(size);
                temp.device(tenx::omp::getDevice()) = mps->get_M_bare().slice(off, ext).contract(ampl, tenx::idx({2}, {0})).reshape(std::array<long, 1>{size});
                t_con.toc();
                ampl = temp;
                // Update the current site
                site = pos;

                // Add to cache
            }
            if(not site) throw std::logic_error(fmt::format("ampl has undefined site: should be {}", tgt_pos));
            if(site.value() != tgt_pos) throw std::logic_error(fmt::format("site ({}) != mps_site ({})", site.value(), tgt_pos));
        }
    }

    void eval(const StateFinite &state, long tgt_pos, std::vector<Amplitude> &cache) {
        if(tgt_pos <= state.get_position<long>())
            eval_from_A(state, tgt_pos, cache);
        else
            eval_from_B(state, tgt_pos, cache);
    }
};

std::vector<double> compute_probability(const StateFinite &state, long tgt_pos, std::vector<Amplitude> &amplitudes, std::vector<Amplitude> &cache) {
    // Here we compute the probability of finding

    auto                t_prob    = tid::tic_scope("probability");
    auto                state_pos = state.get_position<long>();
    auto                state_len = state.get_length<long>();
    auto                tgt_rpos  = state_len - 1 - tgt_pos;
    auto                prob_size = tgt_pos > state_pos ? tgt_rpos + 2 : tgt_pos + 2;
    std::vector<double> probability(static_cast<size_t>(prob_size), 0.0);
    double              probability_sum = 0.0;
    // Figure out which schmidt values to use
    auto                     t_figout = tid::tic_scope("figout");
    Eigen::Tensor<double, 1> schmidt_values;
    if(tgt_pos < state_pos)
        schmidt_values = state.get_mps_site(tgt_pos + 1).get_L().abs(); // A-site
    else if(tgt_pos == state_pos)
        schmidt_values = state.get_mps_site(tgt_pos).get_LC().abs(); // AC-site
    else {
        const auto &mps_left = state.get_mps_site(tgt_pos - 1);
        schmidt_values       = mps_left.isCenter() ? mps_left.get_LC().abs() : mps_left.get_L().abs(); // B-site
    }
    t_figout.toc();

    // Create optional slots for each schmidt value
    auto   t_slots = tid::tic_scope("slots");
    double cutoff  = 1e-14;
    for(long alpha = 0; alpha < schmidt_values.size(); alpha++) {
        auto sq = schmidt_values[alpha] * schmidt_values[alpha];
        for(auto &a : amplitudes) {
            auto n = a.bits.count();
            a.eval(state, tgt_pos, cache); // Evaluate the amplitude
            if(a.ampl.size() != schmidt_values.size()) {
                tools::log->dump_backtrace();
                throw std::logic_error(fmt::format("Mismatching size ampl {} != {}", a.ampl.size(), schmidt_values.size()));
            }
            auto ampl_sqr = std::abs(std::conj(a.ampl[alpha]) * a.ampl[alpha]);
            // Check that the probability would not grow too large, in case we are erroneously considering an amplitude
            // This is important when we work with a small cutoff, where sometimes numerical noise is mistaken for a signal.
            if(ampl_sqr > cutoff and probability_sum + sq <= 1.0 + 1e-8) {
                probability[n] += sq;
                probability_sum += sq;
                if constexpr(settings::debug_numen)
                    tools::log->trace(FMT_STRING("site {:>2} | n {:>2} | bits {} | cutoff {:8.2e} | 1-P {:8.2e} | a({:>4})² {:8.2e} | λ({:>4})² {:8.2e} *"),
                                      tgt_pos, n, a.to_string(tgt_pos <= state_pos), cutoff, 1 - probability_sum, alpha, ampl_sqr, alpha, sq);
                break;

            } else {
                if constexpr(settings::debug_numen)
                    tools::log->trace(FMT_STRING("site {:>2} | n {:>2} | bits {} | cutoff {:8.2e} | 1-P {:8.2e} | a({:>4})² {:8.2e} | λ({:>4})² {:8.2e}"),
                                      tgt_pos, n, a.to_string(tgt_pos <= state_pos), cutoff, 1 - probability_sum, alpha, ampl_sqr, alpha, sq);
            }
        }
    }

    // Sanity check on probabilities
    auto p_sum = std::accumulate(probability.begin(), probability.end(), 0.0);
    if(std::abs(p_sum - 1.0) > 1e-8) {
        tools::log->dump_backtrace();
        tools::log->info("p(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
        throw std::runtime_error(fmt::format("p_sum - 1.0 = {:.8e}", p_sum - 1.0));
    }
    tools::log->trace("p(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
    return probability;
}

std::vector<Amplitude> generate_amplitude_list(const StateFinite &state, long mps_pos) {
    // Generate a list of bit sequences of size prod_{i=0}^pos spin_dim_i.
    // For spin-half this is just 2^pos elements, where pos is the mps position
    // counting from the left.
    // Example: Let mpo_pos == 2, state_pos == 4 and state_len == 8
    //      Then num_bitseqs = prod_{i=0}^{mpo_pos} = spin_dim_0 * spin_dim_1 * spin_dim_2 = 2³ = 8
    //      if mpo_pos == 2 <= state_pos we generate
    //          000, 100, 010, 110, 001, 101, 011, 111
    //      if mpo_pos == 6 > state_pos we generate
    //          000, 100, 010, 110, 001, 101, 011, 111
    // So we generate the same regardless, but later on,
    // these are read in reverse in the second case

    auto t_amp     = tid::tic_scope("amplitude");
    auto state_pos = state.get_position<long>();
    auto state_len = state.get_length<long>();
    auto spinprod  = [](long &acc, const auto &mps) {
        if(mps->spin_dim() != 2) throw std::runtime_error("number_entropies: spin_dim() != 2 is not supported");
        return acc * mps->spin_dim();
    };

    long num_bitseqs;
    if(mps_pos <= state_pos) {
        num_bitseqs = std::accumulate(state.mps_sites.begin(), state.mps_sites.begin() + mps_pos + 1, 1l, spinprod);
    } else {
        auto mps_rpos = state_len - 1 - mps_pos;
        num_bitseqs   = std::accumulate(state.mps_sites.rbegin(), state.mps_sites.rbegin() + mps_rpos + 1, 1l, spinprod);
    }

    std::vector<Amplitude> amplitudes;
    amplitudes.reserve(static_cast<size_t>(num_bitseqs));
    for(long count = 0; count < num_bitseqs; count++) amplitudes.emplace_back(Amplitude{state_len, std::bitset<64>(count), std::nullopt, {}});
    return amplitudes;
}

/*! \brief Compute the number entropies at all cuts of the 1d-chain
 *
 * In fLBIT simulations we have particle-number conservation.
 * The number of particles is invariant to both time evolution of the lbit state,
 * as well as transforming to real-basis with the unitary circuit.
 *
 * The number entropy is defined as
 *
 *      S_N = - Σ_n p(n) ln p(n)
 *
 * where p(n) is the probability of finding n particles in a subsystem l.
 * At a given bond, we calculate p(n) for the left/right subsystem
 * from schmidt values corresponding to nonzero amplitudes A(n)_𝛼:
 *      p(n) =  Σ_[𝛼|A({𝜎})_𝛼 > ϵ] (λ_α)²,
 *
 * where:
 *      * {𝜎} is a particle configuration on the 1d chain such as 11000 or 00101
 *      * n is the number of particles in {𝜎},  Σ_i 𝜎_i.  (e.g. both 11000 and 00101 would give n = 2)
 *      * A({𝜎}) is the amplitude vector identified by the configuration {𝜎}. There are 2^l per bond.
 *      * λ are the schmidt values of the bond at the bipartition. There is only one λ per bond.
 *      * 𝛼 indexes both the schmidt values λ the amplitude vector and A({𝜎})
 *      * ϵ is a small cutoff, typically ~ 1e-14
 *      * subscript _[𝛼|A({𝜎})_𝛼 > ϵ] means that the sum only includes
 *      * 𝛼 for which the corresponding amplitude element is large enough, signaling that
 *        the probability of finding n particles in the subsystem is nonzero.
 *
 *
 * Example: Consider a bipartition in the middle of 8 particles.
 *          In general, an MPS describing such a subsystem has the form
 *
 *            [0]--[1]--[2]--[3]--A({𝜎})_α
 *             |    |    |    |
 *            𝜎0   𝜎1   𝜎2   𝜎3
 *
 *          where 𝜎i = [0,1] sets an lbit off/on on that site,
 *          and 𝛼 indexes both the A({𝜎}) and schmidt values λ.
 *          Observe that:
 *              * Each combination of 𝜎i's such as {𝜎}=01101 yields a unique amplitude vector A(01101)
 *              * We expect A({𝜎}) and A({𝜎'}) with the same n to have nonzero elements at the same α
 *
 *
 */
std::vector<double> tools::finite::measure::number_entropies(const StateFinite &state) {
    if(state.measurements.number_entropies) return state.measurements.number_entropies.value();
    if(state.get_algorithm() != AlgorithmType::fLBIT) {
        // Only fLBIT has particle-number conservation
        throw std::logic_error(
            fmt::format("Called number_entropies(StateFinite) from algorithm [{}]. Only [fLBIT] is allowed.", enum2sv(state.get_algorithm())));
    }

    auto t_num      = tid::tic_scope("number_entropy");
    auto state_copy = state; // Make a local copy so we can move it to the middle without touching the original state
    tools::finite::mps::move_center_point_to_middle(state_copy, state_copy.find_largest_chi());

    auto state_pos       = state_copy.get_position<long>();
    auto state_len       = state_copy.get_length();
    auto von_neumann_sum = [](double sum, const double p) {
        return p > 0 ? sum + p * std::log(p) : sum;
    };

    std::vector<double>                                 number_entropies(state_len + 1, 0.0); // Collects the resulting number entropies
    std::vector<Amplitude>                              cache;
    std::vector<std::pair<size_t, std::vector<double>>> probabilities;
    tools::log->enable_backtrace(200);
    for(const auto &mps : state_copy.mps_sites) {
        auto pos = mps->get_position<long>();
        auto idx = static_cast<size_t>(pos) + 1; // First [0] and last [L+1] number entropy are zero. Then mps[0] generates number entropy idx 1, and so on.
        if(pos > state_pos) break;               // Only compute up to and including AC
        if(mps->get_label() == "B") throw std::logic_error("Expected A/AC site, got B");
        auto amplitudes       = generate_amplitude_list(state_copy, pos);
        auto probability      = compute_probability(state_copy, pos, amplitudes, cache);
        auto number_entropy   = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
        number_entropies[idx] = std::abs(number_entropy);
        cache                 = amplitudes; // Cache the amplitudes for the next step
        if constexpr(settings::debug_numen) probabilities.emplace_back(idx, probability);
    }

    cache.clear();

    for(const auto &mps : iter::reverse(state_copy.mps_sites)) {
        auto pos = mps->get_position<long>();
        auto idx = static_cast<size_t>(pos); // First [0] and last [L+1] number entropy are zero. Then mps[L] generates number entropy idx L, and so on.
        if(pos <= state_pos + 1) break;      // No need to compute at AC again so add +1
        if(mps->get_label() != "B") throw std::logic_error(fmt::format("Expected B site, got {}", mps->get_label()));
        auto amplitudes       = generate_amplitude_list(state_copy, pos);
        auto probability      = compute_probability(state_copy, pos, amplitudes, cache);
        auto number_entropy   = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
        number_entropies[idx] = std::abs(number_entropy);
        cache                 = amplitudes; // Cache the amplitudes for the next step
        if constexpr(settings::debug_numen) probabilities.emplace_back(idx, probability);
    }
    if constexpr(settings::debug_numen)
        for(const auto &[idx, prob] : probabilities) {
            // Sanity check on probabilities
            auto p_sum = std::accumulate(prob.begin(), prob.end(), 0.0);
            tools::log->trace("idx {:>2} | p(n) = {:20.16f} = {:20.16f}", idx, fmt::join(prob, ", "), p_sum);
        }

    tools::log->disable_backtrace();
    state.measurements.number_entropies        = number_entropies;
    state.measurements.number_entropy_midchain = number_entropies.at(state.get_length() / 2);
    state.measurements.number_entropy_current  = number_entropies.at(state.get_position<size_t>() + 1);
    tools::log->debug(FMT_STRING("Number entropies: {:.4f}"), fmt::join(number_entropies, ", "));
    return state.measurements.number_entropies.value();
}

double tools::finite::measure::number_entropy_current(const StateFinite &state) {
    if(state.measurements.number_entropy_current) return state.measurements.number_entropy_current.value();
    if(state.measurements.number_entropies) return state.measurements.number_entropies->at(state.get_position<size_t>() + 1);
    auto entropies = number_entropies(state);
    if(not entropies.empty()) {
        state.measurements.number_entropy_current = entropies.at(state.get_position<size_t>() + 1);
        return state.measurements.number_entropy_current.value();
    } else
        return 0;
}

double tools::finite::measure::number_entropy_midchain(const StateFinite &state) {
    if(state.measurements.number_entropy_midchain) return state.measurements.number_entropy_midchain.value();
    if(state.measurements.number_entropies) return state.measurements.number_entropies->at(state.get_length() / 2);
    auto entropies = number_entropies(state);
    if(not entropies.empty()) {
        state.measurements.number_entropy_midchain = entropies.at(state.get_length() / 2);
        return state.measurements.number_entropy_midchain.value();
    } else
        return 0;
}