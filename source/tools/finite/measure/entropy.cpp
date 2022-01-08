
#include "../measure.h"
#include "debug/exceptions.h"
#include "math/linalg/matrix.h"
#include <bitset>
#include <config/debug.h>
#include <general/iter.h>
#include <math/num.h>
#include <math/tenx.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/mps.h>
#include <utility>

namespace settings {
    inline constexpr bool debug_numen = false;
}

struct Amplitude {
    long                                  state_size;          // System size
    std::bitset<64>                       bits;                // Bits that select spins on each MPS site
    std::optional<long>                   site = std::nullopt; // MPS site (not Schmidt site!)
    Eigen::Tensor<StateFinite::Scalar, 1> ampl;                // Accumulates the MPS tensors
    bool                                  cache_hit = false;   // True if eval was avoided due to cache hit
    Amplitude(long state_size, std::bitset<64> bits, std::optional<long> site, Eigen::Tensor<StateFinite::Scalar, 1> ampl)
        : state_size(state_size), bits(bits), site(site), ampl(std::move(ampl)) {}

    public:
    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_string(const std::bitset<64> &b, long num) const {
        if constexpr(on) {
            if(tools::log->level() > spdlog::level::trace) return {};
            assert(num > 0 and num < 64 and "num should be in range [0,64]");
            return fmt::format(FMT_STRING("{1:>{0}}"), state_size, b.to_string().substr(b.size() - static_cast<size_t>(num)));
        } else
            return {};
    }
    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_rstring(const std::bitset<64> &b, long num) const {
        if constexpr(on) {
            if(tools::log->level() > spdlog::level::trace) return {};
            if(num < 0l or num::cmp_greater_equal(num, b.size())) throw std::logic_error(fmt::format("num should be in range [0,{}]", b.size()));
            auto bs = b.to_string().substr(b.size() - static_cast<size_t>(num));
            return fmt::format(FMT_STRING("{1:<{0}}"), state_size, std::string{bs.rbegin(), bs.rend()});
        } else
            return {};
    }

    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_string() const {
        if(site)
            return to_string<on>(bits, state_size);
        else
            return {};
    }
    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_rstring() const {
        if(site)
            return to_rstring<on>(bits, state_size);
        else
            return {};
    }
    template<bool on = settings::debug_numen>
    [[nodiscard]] std::string to_string(bool reverse) const {
        if(site)
            if(reverse)
                return to_rstring<on>(bits, state_size);
            else
                return to_string<on>(bits, state_size);
        else
            return {};
    }

    [[nodiscard]] bool contains(const Amplitude &other, long tgt_pos) const {
        // Example. Let this amplitude be a and the other is o.
        // Consider a bitset of 8 bits, and we are checking tgt_pos == 3:
        // Then bits 0,1 and 2 need to be identical
        //     a.bits = 00001011
        //     o.bits = 00000011
        // Shift bits to the left edge with <<
        //     a << (8 - tgt_pos) = 01100000
        //     o << (8 - tgt_pos) = 01100000
        // Compare with the ^ operator and check that all are zero with .none()
        //     r = (( a << (8 - tgt_pos) ) ^ (o << (8 - tgt_pos))) = 00000000
        //     return r.none()
        const auto &a   = this->bits;
        const auto &o   = other.bits;
        const auto  pos = static_cast<size_t>(tgt_pos);
        return ((a << (a.size() - pos)) ^ (o << (o.size() - pos))).none();
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
                        cache_hit  = c.site.value() + 1 == tgt_pos;
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
                cache_hit = false;
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
                        cache_hit  = c.site.value() + 1 == tgt_pos;
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
                cache_hit = false;
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

size_t nextGreaterWithSameSetBit(size_t n) {
    // http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    //     size_t t = n | (n - 1); // t gets v's least significant 0 bits set to 1
    //  Next set to 1 the most significant bit to change,
    //  set to 0 the least significant ones, and add the necessary 1 bits.
    //     return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(n) + 1));

    // The __builtin_ctz(v) GNU C compiler intrinsic for x86 CPUs returns the number of trailing zeros. If you are using Microsoft compilers for x86, the
    // intrinsic is _BitScanForward. These both emit a bsf instruction, but equivalents may be available for other architectures. If not, then consider using
    // one of the methods for counting the consecutive zero bits mentioned earlier.
    //
    // Here is another version that tends to be slower because of its division operator, but it does not require counting the trailing zeros.

    if(n == 0ul) return -1ul;
    size_t t = (n | (n - 1)) + 1;
    return t | ((((t & -t) / (n & -n)) >> 1) - 1);
}

// function to find the position of rightmost set bit. Returns -1 if there are no set bits
size_t getFirstSetBitPos(size_t n) { return static_cast<size_t>((std::log2(n & -n) + 1) - 1); }

// function to find the next greater integer
size_t nextGreaterWithOneMoreSetBit(size_t n) {
    // position of rightmost unset bit of n by passing ~n as argument
    size_t pos = getFirstSetBitPos(~n);
    // if n consists of unset bits, then set the rightmost unset bit
    if(pos > -1ul) return (1 << pos) | n;
    // n does not consists of unset bits
    return ((n << 1) + 1);
}

size_t amplitude_next_idx_round_robin(size_t aidx, size_t &nbit, size_t nmax, std::vector<bool> &nflg, std::vector<size_t> &namp,
                                      const std::vector<Amplitude> &amplitudes) {
    if(namp.size() != amplitudes.size()) throw except::logic_error("Size mistmatch: namp {} != amplitudes {}", namp.size(), amplitudes.size());
    if(nflg.size() != nmax + 1) throw except::logic_error("nflg.size() != nmax+1");
    auto t_next = tid::tic_token("next");
    // aidx is the current index in abit
    // nbit is the number of bits at this index, i.e. nbit = abit[aidx].
    if(nbit == namp[aidx]) {
        namp[aidx] = -1ul; // Mark as used
        return aidx;       // The current index hasn't been checked, so return it
    }

    auto next_nbit = [&](size_t n) -> std::pair<size_t, bool> {
        if(nflg.size() != nmax + 1) throw except::logic_error("nflg.size() != nmax+1");
        bool found_all_nbit = std::all_of(nflg.begin(), nflg.end(), [](bool b) -> bool { return b; });
        if(found_all_nbit) return {-1ul, found_all_nbit};
        while(true) {
            n = num::pbc<size_t>(n + 1, nmax + 1); // Advance
            if(not nflg[n]) break;                 // Keep if it's still valid
        }
        return {n, found_all_nbit};
    };

    // The goal now is to find the next index aidx in abit which has nbit+1
    size_t counter                 = 0;
    size_t wcount                  = 0;
    bool   found_all_nbit          = false;
    std::tie(nbit, found_all_nbit) = next_nbit(nbit);

    while(not found_all_nbit) {
        wcount++;
        if(nbit == namp[aidx]) {
            namp[aidx] = -1ul; // Mark as used
            counter    = 0;    // Reset counter when we find a match
                               //            tools::log->info("aidx {} | nbit {} | wcount: {}",aidx, nbit, wcount);
            return aidx;       // Return aidx. This index points to the next amplitude that has nbit increased by one.
        }

        if(namp[aidx] == -1ul and nbit == std::bitset<64>(aidx).count()) {
            // aidx has the correct nbit but this particular aidx is already taken.
            // We can fast-forward to the next aidx with the same nbit and try that
            auto aidx_new = nextGreaterWithSameSetBit(aidx);
            if(aidx_new < namp.size()) {
                counter += aidx_new - aidx;
                aidx = aidx_new;
                continue; // Try this aidx instead
            }
        }
        aidx = num::pbc<size_t>(aidx + 1, namp.size()); // Cycle through the whole list of amplitudes

        if(counter++ >= namp.size()) {
            // We haven't been able to find a matching index, which means that all indices with the current nbit are taken.
            // Simply go to the next nbit and reset the counter
            counter                        = 0;
            nflg[nbit]                     = true; // Flag this nbit as used up
            auto nbit_old                  = nbit;
            std::tie(nbit, found_all_nbit) = next_nbit(nbit);
            if(not found_all_nbit and nbit > nbit_old) {
                // Fastforward to the next greater with one nbit increased appropriately
                auto aidx_new = aidx;
                for(size_t i = 0; i < nbit - nbit_old; i++) aidx_new = nextGreaterWithOneMoreSetBit(aidx_new);
                if(aidx_new < namp.size()) aidx = aidx_new; // Take it if it's in range. Otherwise ignore.
            }
        }
    }
    return -1ul;
}

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
    auto                t_slots            = tid::tic_scope("slots");
    double              amplitude_cutoff   = 1e-12;
    double              probability_cutoff = 1e-32;
    std::vector<size_t> namp; // Number of bits in each amplitude
    namp.reserve(amplitudes.size());
    for(const auto &a : amplitudes) namp.emplace_back(a.bits.count());

    auto              nmax = static_cast<size_t>(tgt_pos) + 1ul; // Maximum number of n
    std::vector<bool> nflg(nmax + 1, false);                     // Flags for each n to tell if they have all been found in amplitudes

    Eigen::MatrixXd ampacc_sq_matrix = Eigen::MatrixXd::Zero(schmidt_values.size(), tgt_pos + 2);
    Eigen::VectorXi schmidt_taken    = Eigen::VectorXi::Zero(schmidt_values.size());
    Eigen::VectorXd schmidt_squared  = tenx::VectorMap(schmidt_values).cwiseAbs2();
    auto            idx              = 0ul; // Start amplitude index
    auto            n                = 0ul; // Number of bits in the current amplitude
    auto            min_alpha        = 0l;
    auto            max_alpha        = schmidt_values.size();

    while(true) {
        idx = amplitude_next_idx_round_robin(idx, n, nmax, nflg, namp, amplitudes);
        if(idx == -1ul) // Could not find next idx. Probably all have been checked.
            break;
        auto &a  = amplitudes[idx];
        long  nl = static_cast<long>(a.bits.count());
        if constexpr(settings::debug_numen)
            if(n != nl) throw except::logic_error("Wrong bit number!");
        // Evaluate the amplitude vector
        a.eval(state, tgt_pos, cache);
        if(a.ampl.size() != schmidt_values.size()) {
            tools::log->dump_backtrace();
            throw std::logic_error(fmt::format("Mismatching size ampl {} != {}", a.ampl.size(), schmidt_values.size()));
        }
        // Add amplitudes to the nth column
        auto avec = tenx::VectorMap(a.ampl);
        ampacc_sq_matrix.col(nl) += avec.conjugate().cwiseProduct(avec).cwiseAbs();
        //        tools::log->info("ampacc_sq: \n{}\n", linalg::matrix::to_string(ampacc_sq_matrix,8));

        // Check if any amplitude element gives the signal to add probability
        for(long alpha = min_alpha; alpha < max_alpha; alpha++) {
            if(schmidt_taken(alpha) == 1) continue;
            auto asq = ampacc_sq_matrix(alpha, nl); // The a²[alpha] value tells us to pick the corresponding λ²[alpha] when nonzero.
            auto ssq = schmidt_squared[alpha];      // The value λ² is added to probability if the amplitude is greater than cutoff.
            // Check that the probability would not grow too large, in case we are erroneously considering an amplitude
            // This is important when we work with a small cutoff, where sometimes numerical noise is mistaken for a signal.
            bool accept = asq > amplitude_cutoff and probability_sum + ssq <= 1.0 + 1e-8;
            if constexpr(settings::debug_numen) {
                char accept_ch = accept ? '*' : ' ';
                char cacheh_ch = a.cache_hit ? 'o' : ' ';
                tools::log->trace(
                    FMT_STRING("site {:>2} | n {:>2} | bits {} | idx {:>5} | cutoff {:8.2e} | 1-P {:8.2e} | a({:>4})² {:8.2e} | λ({:>4})² {:8.2e} [{}|{}]"),
                    tgt_pos, n, a.to_string(tgt_pos <= state_pos), idx, amplitude_cutoff, 1 - probability_sum, alpha, asq, alpha, ssq, accept_ch, cacheh_ch);
            }

            if(accept) {
                probability[n] += ssq;
                probability_sum += ssq;
                schmidt_taken(alpha) = 1;
            }
        }
        if(schmidt_taken.isOnes()) break;                                                    // All schmidt values squared have been added to probability
        while(schmidt_taken(min_alpha) == 1) min_alpha++;                                    // Advance min_alpha to skip first taken schmidt values
        while(schmidt_taken(max_alpha - 1) == 1 and max_alpha >= min_alpha + 1) max_alpha--; // Decrease max_alpha to skip the last taken schmidt values
    }

    // Sanity check on probabilities
    auto p_sum = std::accumulate(probability.begin(), probability.end(), 0.0);
    if(std::abs(p_sum - 1.0) > 1e-4) {
        tools::log->dump_backtrace();
        tools::log->info("p(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
        throw except::runtime_error("p_sum - 1.0 = {:.8e}", p_sum - 1.0);
    }
    if(std::abs(p_sum - 1.0) > 1e-8) {
        tools::log->dump_backtrace();
        tools::log->warn("p(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
        tools::log->warn("p_sum - 1.0 = {:.8e}", p_sum - 1.0);
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
    // So we generate the same regardless, but later on these are read in reverse in the second case

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
    for(long count = 0; count < num_bitseqs; count++)
        amplitudes.emplace_back(Amplitude{state_len, std::bitset<64>(static_cast<unsigned long long int>(count)), std::nullopt, {}});

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
 *        𝛼 for which the corresponding amplitude element is large enough, signaling that
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

    if constexpr(settings::debug_numen) {
        for(const auto &[idx, prob] : probabilities) {
            // Sanity check on probabilities
            auto p_sum = std::accumulate(prob.begin(), prob.end(), 0.0);
            tools::log->trace("idx {:>2} | p(n) = {:20.16f} = {:20.16f}", idx, fmt::join(prob, ", "), p_sum);
        }
    }

    tools::log->disable_backtrace();
    state.measurements.number_entropies        = number_entropies;
    state.measurements.number_entropy_midchain = number_entropies.at(state.get_length<size_t>() / 2);
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
