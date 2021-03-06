
#include "../measure.h"
#include <bitset>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_extra.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/mps.h>

struct Amplitude {
    std::bitset<64>                              bits;                // Bits that select spins on each MPS site
    std::optional<long>                          site = std::nullopt; // MPS site (not Schmidt site!)
    Eigen::Tensor<class_state_finite::Scalar, 1> ampl;                // Accumulates the MPS tensors

    [[nodiscard]] std::string to_string(const std::bitset<64> &b, long num) const {
        std::string s;
        if(num > 64) throw std::range_error("Maximum bitset length is 64");
        for(size_t i = 0; i <= static_cast<size_t>(num); i++) s.append(std::to_string(static_cast<int>(b[i])));
        return s;
    }
    [[nodiscard]] std::string to_string() const {
        if(site)
            return to_string(bits, site.value());
        else
            return std::string();
    }
    void eval_from_A(const class_state_finite &state, long tgt_pos, const std::vector<Amplitude> &cache = {}) {
        // Start by calculating the mps sites that should be included in the amplitude
        if(not site or site.value() < tgt_pos) {
            // There are missing mps in the amplitude. Let's add them
            // First, check if it is already available in the amplitudes database
            // We can look it up by transforming the bits between 0 to bit_site to an integer
            // and use that as the index on amplitudes
            long state_pos = state.get_position<long>();
            if(tgt_pos > state_pos) throw std::logic_error(fmt::format("eval_from_A: expected mps_site ({}) <= state_pos ({})", tgt_pos, state_pos));

            if(tgt_pos > 0 and not cache.empty()) {
                // There is a chance to continue building an existing amplitude
                std::bitset<64> bits_index = 0;
                for(size_t i = 0; i < static_cast<size_t>(tgt_pos); i++) bits_index[i] = bits[i];
                if(cache.size() > bits_index.to_ulong()){
                    const auto &c = cache.at(bits_index.to_ulong()); // A cache item
                    if(c.site and c.site.value() <= tgt_pos) {
                        // Cache hit! No need to compute the amplitude from scratch
                        tools::log->trace("Cached   site {} | bits [{}] at idx {}", c.site.value(), c.to_string(), bits_index.to_ulong());
                        site = c.site;
                        ampl = c.ampl;
                    }
                }
            }

            if(ampl.size() == 0) { // Initialize
                ampl.resize(1);
                ampl.setConstant(1.0);
            }
            // Contract the missing mps up to, but not including, the last mps at mps_site
            for(const auto &mps : state.mps_sites) {
                long pos = mps->get_position<long>();
                if(site and site.value() >= pos) continue; // Fast-forward to the missing sites
                if(tgt_pos < pos) break;                   // Contract up to the mps at tgt_pos
                if(ampl.size() != mps->get_chiL())
                    throw std::runtime_error(fmt::format("eval() failed for site {}: "
                                                         "mismatch in ampl({}) with size = {} and mps({}) with chiL = {} | bits {}",
                                                         tgt_pos, site.value(), ampl.size(), pos, mps->get_chiL(), bits.to_string()));
                long                size = mps->get_chiR();
                std::array<long, 3> off  = {bits[static_cast<size_t>(pos)], 0, 0};
                std::array<long, 3> ext  = {1, mps->get_chiL(), mps->get_chiR()};

                // ampl never has a trailing Lambda, which means that we must SVD lambda out of B sites
                Eigen::Tensor<class_state_finite::Scalar, 1> temp =
                    ampl.contract(mps->get_M_bare().slice(off, ext), Textra::idx({0}, {1})).reshape(std::array<long, 1>{size});
                ampl = temp;

                // Update the current site
                site = pos;
            }
            if(not site) throw std::logic_error(fmt::format("ampl has undefined site: should be {}", tgt_pos));
            if(site.value() != tgt_pos) throw std::logic_error(fmt::format("site ({}) != mps_site ({})", site.value(), tgt_pos));
        }
    }

    void eval_from_B(const class_state_finite &state, long tgt_pos, const std::vector<Amplitude> &cache) {
        // Start by calculating the mps sites that should be included in the amplitude
        // Remember that in eval_from_B we calculate the amplitude starting from the right-end of the chain
        if(not site or site.value() > tgt_pos) {
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
                std::bitset<64> bits_index = 0;
                for(size_t i = 0; i < static_cast<size_t>(mps_rsite); i++) bits_index[i] = bits[i];
                if(cache.size() > bits_index.to_ulong()){
                    const auto &c = cache.at(bits_index.to_ulong()); // A cache item
                    if(c.site and c.site.value() > tgt_pos) {
                        // Cache hit! No need to compute the amplitude from scratch
                        tools::log->trace("Cached   site {} | bits [{}] at idx {}", c.site.value(), c.to_string(), bits_index.to_ulong());
                        site = c.site;
                        ampl = c.ampl;
                    }
                }
            }

            if(ampl.size() == 0) { // Initialize
                ampl.resize(1);
                ampl.setConstant(1.0);
            }
            // Contract the missing mps
            for(const auto &mps : iter::reverse(state.mps_sites)) {
                long pos = mps->get_position<long>();
                if(site and site.value() <= pos) continue; // Fast-forward to the missing sites
                if(tgt_pos > pos) break;                   // Contract up to the mps at mps_pos
                if(ampl.size() != mps->get_chiR())
                    throw std::runtime_error(fmt::format("eval() failed for site {}: "
                                                         "mismatch in ampl({}) with size = {} and mps({}) with chiR = {} | bits {}",
                                                         tgt_pos, site.value(), ampl.size(), pos, mps->get_chiL(), bits.to_string()));
                long                mps_rpos = state_len - 1 - mps->get_position<long>();
                long                size     = mps->get_chiL();
                std::array<long, 3> off      = {bits[static_cast<size_t>(mps_rpos)], 0, 0};
                std::array<long, 3> ext      = {1, mps->get_chiL(), mps->get_chiR()};

                // ampl never has a trailing Lambda, which means that we must SVD lambda out of B sites
                Eigen::Tensor<class_state_finite::Scalar, 1> temp =
                    mps->get_M_bare().slice(off, ext).contract(ampl, Textra::idx({2}, {0})).reshape(std::array<long, 1>{size});
                ampl = temp;
                // Update the current site
                site = pos;
            }
            if(not site) throw std::logic_error(fmt::format("ampl has undefined site: should be {}", tgt_pos));
            if(site.value() != tgt_pos) throw std::logic_error(fmt::format("site ({}) != mps_site ({})", site.value(), tgt_pos));
        }
    }

    void eval(const class_state_finite &state, long tgt_pos, const std::vector<Amplitude> &database = {}) {
        if(tgt_pos <= state.get_position<long>())
            eval_from_A(state, tgt_pos, database);
        else
            eval_from_B(state, tgt_pos, database);
    }
};

std::vector<double> compute_probability(const class_state_finite &state, long tgt_pos, std::vector<Amplitude> &amplitudes,
                                        const std::vector<Amplitude> &cache = {}) {
    auto                state_pos = state.get_position<long>();
    auto                state_len = state.get_length<long>();
    auto                tgt_rpos  = state_len - 1 - tgt_pos;
    auto                prob_size = tgt_pos > state_pos ? tgt_rpos + 2 : tgt_pos + 2;
    std::vector<double> probability(static_cast<size_t>(prob_size), 0.0);

    // Figure out which schmidt values to use
    Eigen::Tensor<double, 1> schmidt_values;
    if(tgt_pos < state_pos)
        schmidt_values = state.get_mps_site(tgt_pos + 1).get_L().abs(); // A-site
    else if(tgt_pos == state_pos)
        schmidt_values = state.get_mps_site(tgt_pos).get_LC().abs(); // AC-site
    else {
        const auto &mps_left = state.get_mps_site(tgt_pos - 1);
        schmidt_values       = mps_left.isCenter() ? mps_left.get_LC().abs() : mps_left.get_L().abs(); // B-site
    }

    // Create optional slots for each schmidt value
    std::vector<double> schmidt_vector(static_cast<size_t>(schmidt_values.size()));
    for(auto &&[i, s] : iter::enumerate(schmidt_vector)) s = schmidt_values[static_cast<long>(i)];
    double cutoff = 1.0 / std::pow(state_len, 2);
    for(const auto &[i, s] : iter::enumerate(schmidt_vector)) {
        for(auto &a : amplitudes) {
            // Evaluate the amplitude
            a.eval(state, tgt_pos, cache);
            if(a.ampl.size() != static_cast<long>(schmidt_vector.size())) {
                tools::log->dump_backtrace();
                throw std::runtime_error(fmt::format("Mismatching size ampl {} != {}", a.ampl.size(), schmidt_vector.size()));
            }
            auto ampl_sqr = std::abs(std::conj(a.ampl[static_cast<long>(i)]) * a.ampl[static_cast<long>(i)]);
            if(ampl_sqr > cutoff) {
                auto n  = a.bits.count();
                auto sq = s * s;
                probability[n] += sq;
                tools::log->trace("site {} | n {} | bits {}  | i {} | c {:22.20f} | ampl_sqr = {:22.20f} *", tgt_pos, a.bits.count(), a.to_string(), i, cutoff,
                                  ampl_sqr);
                break;
            }
//            else {
//                tools::log->trace("site {} | n {} | bits {}  | i {} | c {:22.20f} | ampl_sqr = {:22.20f}", tgt_pos, a.bits.count(), a.to_string(), i, cutoff,
//                                  ampl_sqr);
//            }
        }
    }

    // Sanity check on probabilities

    auto p_sum = std::accumulate(probability.begin(), probability.end(), 0.0);
    if(std::abs(p_sum - 1.0) > 1e-8) {
        tools::log->dump_backtrace();
        tools::log->info("p(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
        throw std::runtime_error(fmt::format("p_sum - 1.0 = {:.8e}", p_sum - 1.0));
    }
    //    tools::log->info("p(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
    return probability;
}

std::vector<Amplitude> generate_amplitude_list(const class_state_finite &state, long mps_pos) {
    auto state_pos = state.get_position<long>();
    auto state_len = state.get_length<long>();
    auto spinprod  = [](long &acc, const auto &mps) {
        if(mps->spin_dim() != 2) throw std::runtime_error("number_entropies: spin_dim() != 2 is not supported");
        return acc * mps->spin_dim();
    };
    // Generate a list of bit sequences of size prod_{i=0}^l spin_dim_i.
    // For spin-half this is just 2^l elements, where l is the mps position
    long num_bitseqs;
    if(mps_pos <= state_pos) {
        num_bitseqs = std::accumulate(state.mps_sites.begin(), state.mps_sites.begin() + mps_pos + 1, 1l, spinprod);
    } else {
        auto mps_rpos = state_len - 1 - mps_pos;
        num_bitseqs   = std::accumulate(state.mps_sites.rbegin(), state.mps_sites.rbegin() + mps_rpos + 1, 1l, spinprod);
    }

    std::vector<Amplitude> amplitudes;
    amplitudes.reserve(static_cast<size_t>(num_bitseqs));
    for(long count = 0; count < num_bitseqs; count++)
        amplitudes.emplace_back(Amplitude{static_cast<unsigned long long int>(count), std::nullopt, {}});
    return amplitudes;
}

std::vector<double> tools::finite::measure::number_entropies(const class_state_finite &state) {
    if(state.measurements.number_entropies) return state.measurements.number_entropies.value();
    if(state.get_algorithm() != AlgorithmType::fLBIT) {
        // Only fLBIT has particle-number conservation
        state.measurements.number_entropy_current = {};
        return {};
    }

    if(state.get_algorithm() != AlgorithmType::fLBIT) return {}; // Only fLBIT has particle-number conservation
    auto state_copy = state; // Make a local copy so we can move it to the middle without touching the original state
    tools::finite::mps::move_center_point_to_middle(state_copy, state_copy.find_largest_chi());

    auto t_num = tools::common::profile::get_default_prof()["t_num"]->tic_token();
    auto state_pos       = state_copy.get_position<long>();
    auto state_len       = state_copy.get_length();
    auto von_neumann_sum = [](double sum, const double p) {
        return p > 0 ? sum + p * std::log(p) : sum;
    };

    std::vector<double>    number_entropies(state_len + 1, 0.0); // Collects the resulting number entropies
    std::vector<Amplitude> cache;
    tools::log->enable_backtrace(200);
    for(const auto &mps : state_copy.mps_sites) {
        auto pos = mps->get_position<long>();
        if(pos > state_pos) break;
        if(mps->get_label() == "B") throw std::logic_error("Expected A/AC site, got B");
        auto amplitudes                                   = generate_amplitude_list(state_copy, pos);
        auto probability                                  = compute_probability(state_copy, pos, amplitudes, cache);
        auto number_entropy                               = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
        number_entropies[mps->get_position<size_t>() + 1] = std::abs(number_entropy);
        cache                                             = amplitudes; // Cache the amplitudes for the next step
    }

    cache.clear();

    for(const auto &mps : iter::reverse(state_copy.mps_sites)) {
        auto pos = mps->get_position<long>();
        if(pos <= state_pos+1) break; // No need to compute at LC again
        if(mps->get_label() != "B") throw std::logic_error(fmt::format("Expected B site, got {}", mps->get_label()));
        auto amplitudes                               = generate_amplitude_list(state_copy, pos);
        auto probability                              = compute_probability(state_copy, pos, amplitudes, cache);
        auto number_entropy                           = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
        number_entropies[mps->get_position<size_t>()] = std::abs(number_entropy);
        cache                                         = amplitudes; // Cache the amplitudes for the next step
    }

    tools::log->disable_backtrace();
    state.measurements.number_entropies = number_entropies;
    tools::log->info("Number entropies: {}", number_entropies);
    return state.measurements.number_entropies.value();
}

double tools::finite::measure::number_entropy_current(const class_state_finite &state) {
    if(state.measurements.number_entropy_current) return state.measurements.number_entropy_current.value();
    if(state.get_algorithm() != AlgorithmType::fLBIT) {
        // Only fLBIT has particle-number conservation
        state.measurements.number_entropy_current = 0;
        return 0;
    }
    auto pos = state.get_position<long>();
    if(state.measurements.number_entropies){
        if(state.measurements.number_entropies->size() != state.get_length<size_t>()+1)
            throw std::runtime_error(fmt::format("expected number_entropies.size() == lenght+1. Got: ({})", state.measurements.number_entropies->size()));
        return state.measurements.number_entropies->at(static_cast<size_t>(pos + 1));
    }


    auto state_copy = state; // Make a local copy so we can move it to the middle without touching the original state
    tools::finite::mps::move_center_point_to_middle(state_copy, state.find_largest_chi());
    auto von_neumann_sum = [](double sum, const double p) {
      return p > 0 ? sum + p * std::log(p) : sum;
    };
    auto amplitudes                                   = generate_amplitude_list(state_copy, pos);
    auto probability                                  = compute_probability(state_copy, pos, amplitudes);
    state.measurements.number_entropy_current        = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
    return state.measurements.number_entropy_current.value();
}

double tools::finite::measure::number_entropy_midchain(const class_state_finite &state) {
    if(state.measurements.number_entropy_midchain) return state.measurements.number_entropy_midchain.value();
    if(state.measurements.number_entropies) return state.measurements.number_entropies->at(state.get_length()/2);
    if(state.get_algorithm() != AlgorithmType::fLBIT) {
        // Only fLBIT has particle-number conservation
        state.measurements.number_entropy_current = 0;
        return 0;
    }

    auto von_neumann_sum = [](double sum, const double p) {
      return p > 0 ? sum + p * std::log(p) : sum;
    };

    auto pos = state.get_length<long>()/2;
    auto amplitudes                                   = generate_amplitude_list(state, pos);
    auto probability                                  = compute_probability(state, pos, amplitudes);
    state.measurements.number_entropy_midchain        = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
    return state.measurements.number_entropy_midchain.value();
}