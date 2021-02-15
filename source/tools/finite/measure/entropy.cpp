
#include "../measure.h"
#include <bitset>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <math/stat.h>
#include <math/svd.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/mps.h>

double number_entropy(const class_state_finite &state, size_t site, const Eigen::Tensor<class_state_finite::Scalar, 1> &schmidt_values) {
    struct AmplMeta {
        size_t                                       n; // Number of spin ups
        double                                       p; // probability of finding n spin ups P(n) = sum_i {such that n = i } |schmidt_i|^2
        Eigen::Tensor<class_state_finite::Scalar, 1> ampl;
    };

    // Calculate the number of combinations
    long num_comb = 1;
    for(const auto &mps : state.mps_sites) {
        num_comb *= mps->spin_dim();
        tools::log->info("{}[{}] dims {}\n{}", mps->get_label(), mps->get_position(), mps->dimensions(), linalg::tensor::to_string(mps->get_M(), 3, 2));
        if(site == mps->get_position()) break;
    }
    auto combinations = num::range<std::bitset<32>>(0, num_comb);
    //    std::vector<Eigen::Tensor<class_state_finite::Scalar,1>> amplitudes;
    std::vector<AmplMeta> amplitudes;
    //    Eigen::Tensor<double,1> probabilities(static_cast<long>(site + 2));
    std::vector<double> probabilities(site + 2, 0);

    // We will need the squared schmidt values to compute probabilities
    std::vector<double> schmidt_values_sq;
    schmidt_values_sq.reserve(static_cast<unsigned long>(schmidt_values.size()));
    for(long i = 0; i < schmidt_values.size(); i++) schmidt_values_sq.emplace_back(std::abs(schmidt_values[i] * schmidt_values[i]));

    for(const auto &c : combinations) {
        Eigen::Tensor<class_state_finite::Scalar, 2> temp;
        Eigen::Tensor<class_state_finite::Scalar, 2> chain(1, 1);
        chain.setConstant(1.0);
        for(auto &mps : state.mps_sites) {
            //            tools::log->info("Adding {}[{}] | dims {}: \n{}",mps->get_label(), mps->get_position(), mps->dimensions(),
            //            linalg::tensor::to_string(mps->get_M()));
            long                dim0 = 1; // mps->spin_dim();
            long                dimR = mps->get_chiR();
            long                dimL = chain.dimension(0);
            std::array<long, 3> off  = {c[mps->get_position()], 0, 0};
            std::array<long, 3> ext  = {1, mps->get_chiL(), mps->get_chiR()};
            //            tools::log->info("Slice {} \n{}", static_cast<int>(c[mps->get_position()]), linalg::tensor::to_string(mps->get_M().slice(off,ext)));
            temp = chain.contract(mps->get_M().slice(off, ext), Textra::idx({1}, {1})).reshape(std::array<long, 2>{dimL * dim0, dimR});

            chain = temp;
            if(site == mps->get_position()) break;
        }
        // Check if any amplitude in chain is non-zero

        amplitudes.emplace_back(AmplMeta{c.count(), 0.0, chain.reshape(std::array<long, 1>{chain.size()})});
        tools::log->info("Bitset {}: n = {}  ampl = {}", c.to_string(), c.count(), linalg::tensor::to_string(chain, 3, 2));

        double prob = 0;
        for(const auto &[i, s] : iter::enumerate(schmidt_values_sq)) {
            if(std::abs(amplitudes.back().ampl[static_cast<long>(i)]) > 1e-8) prob += s;
        }

        probabilities[c.count()] = prob;
    }

    auto lambda = [](double a, double b) {
        return b > 0.0 ? a + b * std::log(b) : a;
    };
    double number_entropy = -std::accumulate(probabilities.begin(), probabilities.end(), 0.0, lambda);
    //    Eigen::Tensor<double,0> number_entropy =  -probabilities.contract(probabilities.log().eval(), Textra::idx({0}, {0}));

    tools::log->info("L²   = {:.6f}", fmt::join(schmidt_values_sq, ", "));
    tools::log->info("P(n) = {:.6f}", fmt::join(probabilities, ", "));
    tools::log->info("S_N  = {:.6f}", number_entropy);
    return number_entropy;
}

struct Amplitude {
    std::bitset<64>                              bits;                // Bits that select spins on each MPS site
    std::optional<long>                          site = std::nullopt; // MPS site (not Schmidt site!)
    Eigen::Tensor<class_state_finite::Scalar, 1> ampl;                // Accumulates the MPS tensors
    std::optional<double>                        prob = std::nullopt;

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

    size_t find_minimum(const std::vector<double> & diff, const std::vector<double> &ampl_abs){
        // Find the greatest difference in log. This is where the amplitudes drop off into noise
        auto min_logdiff_val = diff[1];
        auto min_logdiff_idx = 1ul;
        for(const auto & [i,a] : iter::enumerate(diff) ){
            if(i <= 1) continue;
            if(a < min_logdiff_val){
                min_logdiff_val = a;
                min_logdiff_idx = i;
            }
            if(ampl_abs[i] < std::numeric_limits<double>::epsilon()) break;
        }
        return min_logdiff_idx;
    }


    std::optional<double> find_cutoff(){
        if(ampl.size() == 0) return std::nullopt;

        // Convert to double
        std::vector<double> ampl_abs(static_cast<size_t>(ampl.size()));
        for(auto && [i,a] : iter::enumerate(ampl_abs)) a = std::abs(ampl[static_cast<long>(i)]);

        // Sort the amplitudes in descending order
        std::sort(ampl_abs.begin(), ampl_abs.end(), std::greater());

        // Find the max element to make sure the amplitudes are large enough
        if(ampl_abs[0] < 1e-8) return std::nullopt; // Not worth dealing with very small amplitudes
        if(ampl_abs.size() == 1){
            return ampl_abs[0];
        }else if(ampl_abs.size() == 2) {
            if(ampl_abs[1] < 1e-8) return ampl_abs[0];
            else return ampl_abs[1];
        }
        else{
            // Take the logarithms of the amplitudes to make it easier to find the drop-off
            std::vector<double> ampl_log;
            for(const auto & a : ampl_abs ){
                if(a == 0) break;
                // Make sure the smallest log is down to epsilon so we don't get a huge
                // step down to log(~0) = -inf which would dominate when trying to find a drop-off.
                ampl_log.emplace_back(
                    std::log10(std::max(a,std::numeric_limits<double>::epsilon())));
            }
            if(ampl_log.empty()){
//                tools::log->trace("cutoff std::nullopt: ampl_log is empty");
                return std::nullopt;
            }

            // Find the greatest difference in log. This is where the amplitudes drop off into noise
            std::vector<double> ampl_logdiff(ampl_log.size());
            std::adjacent_difference(ampl_log.begin(), ampl_log.end(), ampl_logdiff.begin());
            auto min_logdiff_idx = find_minimum(ampl_logdiff, ampl_abs);
            auto min_logdiff_val = ampl_logdiff[min_logdiff_idx];

//            for(const auto & [i,a] : iter::enumerate(ampl_log)) tools::log->trace("ampl {} = {:22.20f} | log {:8.4f} | diff {:8.4f} | min idx {} | min val {}", i, ampl_abs[i], a, ampl_logdiff[i], min_logdiff_idx,min_logdiff_val);
            if(min_logdiff_idx == 0) {
//                tools::log->trace("cutoff std::nullopt: min_logdiff_idx == 0");
                return std::nullopt; // Max amplitude must be smaller than epsilon
            }
            if(min_logdiff_val > -3) {
//                tools::log->trace("cutoff std::nullopt: min_logdiff_val > -3");
                return std::nullopt; // We want a drop at least x1000
            }
//            {
//                // Find the second greatest difference in log. This is to make sure the dropoff is big enough.
//                ampl_logdiff[min_logdiff_idx] = 0.0;
//                auto sec_logdiff_idx = find_minimum(ampl_logdiff, ampl_abs);
//                auto sec_logdiff_val =  ampl_logdiff[sec_logdiff_idx];
//                // Now make sure the dropoff is much bigger than other dropoffs
//                if(min_logdiff_idx != sec_logdiff_idx and std::abs(sec_logdiff_val - min_logdiff_val) < 2){
//                    tools::log->trace("cutoff std::nullopt: distance to second biggest drop is too small: min {} | sec {}",sec_logdiff_val,min_logdiff_val );
//                    return std::nullopt; //Multiple dropoffs of similar size
//                }
//            }

            // The cuttoff can be set to the middle of the cliff
            return 0.5* (ampl_abs.at(min_logdiff_idx-1) + ampl_abs.at(min_logdiff_idx));
        }



    }


    void eval_prob(const Eigen::Tensor<class_state_finite::Scalar, 1> &schmidt_values) {
        if(ampl.size() != schmidt_values.size()) {
            tools::log->dump_backtrace();
            throw std::range_error(fmt::format("ampl size ({}) != schmidt values size ({})", ampl.size(), schmidt_values.size()));
        }
        // Update the probability of this bit sequence

        // The squared amplitudes add up to 1 + noise. We collect schmidt values corresponding to non-zero amplitudes while adding up their squares.
        // Onec the squared amplitudes >= we can stop collecting schmidt values.
        // It is important to collect them in order of decreasing amplitude so we get the most significant elements for sure.

//        Eigen::Tensor<double, 1> ampl_abs = ampl.abs();
        auto cutoff = find_cutoff();
        if(not cutoff) return;
        std::vector<double> ampl_abs(static_cast<size_t>(ampl.size()));
        for(auto && [i,a] : iter::enumerate(ampl_abs)) a = std::abs(ampl[static_cast<long>(i)]);
        double schmidt_sq_sum = 0;
        for(size_t i = 0; i < ampl_abs.size(); i++) {
            auto max_itr = std::max_element(ampl_abs.begin(), ampl_abs.end());
            auto max_ldx = std::distance(ampl_abs.begin(), max_itr);
            auto max_idx = static_cast<size_t>(max_ldx);
            auto max_val = *max_itr;
            auto max_sqr = max_val * max_val;
            auto schmidt_val = std::abs(schmidt_values[max_ldx]);
            auto schmidt_sqr = schmidt_val * schmidt_val;
            bool   accept      = max_val >= cutoff.value();
            tools::log->trace(
                "ampl {:2} = {:22.20f} | ampl² = {:22.20f} | cutoff {:22.20f} | schm = {:22.20f} | schm² = {:22.20f} {}",
                max_idx, max_val, max_sqr, cutoff.value(), schmidt_val, schmidt_sqr, (accept ? "*" : " "));
            if(accept){
                schmidt_sq_sum += schmidt_sqr;
            }
            ampl_abs[max_idx] = 0;
        }

        if(schmidt_sq_sum > 0.0) prob = schmidt_sq_sum; // If nothing was added we don't know the probability (there was no amplitude)
    }

    void eval_from_A(const class_state_finite &state, long mps_site, const std::vector<Amplitude> &database) {
        // Start by calculating the mps sites that should be included in the amplitude
        if(not site or site.value() < mps_site) {


            // There are missing mps in the amplitude. Let's add them
            // First, check if it is already available in the amplitudes database
            // We can look it up by transforming the bits between 0 to bit_site to an integer
            // and use that as the index on amplitudes
            long state_len = state.get_length<long>();
            long state_pos = state.get_position<long>();
            if(mps_site > state_pos ) throw std::logic_error(fmt::format("eval_from_A: expected mps_site ({}) <= state_pos ({})", mps_site,state_pos));

            if(mps_site > 0) {
                // There is a chance to continue building an existing amplitude
                std::bitset<64> bits_index = 0;
                for(size_t i = 0; i < static_cast<size_t>(mps_site); i++) bits_index[i] = bits[i];

                const auto &c = database.at(bits_index.to_ulong()); // A cache item
                if(c.site and c.site.value() <= mps_site) {
                    // Cache hit! No need to compute the amplitude from scratch
                    tools::log->trace("Cached   site {} | bits [{}] at idx {}", c.site.value(), c.to_string(), bits_index.to_ulong());
                    site = c.site;
                    ampl = c.ampl;
                }
            }

            if(ampl.size() == 0) { // Initialize
                ampl.resize(1);
                ampl.setConstant(1.0);
            }
            // Contract the missing mps up to, but not including, the last mps at mps_site
            for(const auto &mps : state.mps_sites) {
                long mps_pos = mps->get_position<long>();
                if(site and site.value() >= mps_pos) continue; // Fast-forward to the missing sites
                if(mps_site < mps_pos) break;                  // Contract up to the mps at mps_site
                if(ampl.size() != mps->get_chiL())
                    throw std::runtime_error(fmt::format("eval() failed for site {}: "
                                                         "mismatch in ampl({}) with size = {} and mps({}) with chiL = {} | bits {}",
                                                         mps_site, site.value(), ampl.size(), mps_pos, mps->get_chiL(), bits.to_string()));
                long                size = mps->get_chiR();
                std::array<long, 3> off  = {bits[static_cast<size_t>(mps_pos)], 0, 0};
                std::array<long, 3> ext  = {1, mps->get_chiL(), mps->get_chiR()};

                // ampl never has a trailing Lambda, which means that we must SVD lambda out of B sites
                Eigen::Tensor<class_state_finite::Scalar, 1> temp =
                    ampl.contract(mps->get_M_bare().slice(off, ext), Textra::idx({0}, {1})).reshape(std::array<long, 1>{size});
                ampl = temp;

                // Update the current site
                site = mps_pos;
            }
            if(not site) throw std::logic_error(fmt::format("ampl has undefined site: should be {}", mps_site));
            if(site.value() != mps_site) throw std::logic_error(fmt::format("site ({}) != mps_site ({})",site.value(), mps_site ));
            // Update the probability of this bit sequence

            // Now we want "ampl" to end in "... Lambda * Gammma" and no trailing Lambda so that we can study the amplitudes without
            // problems of double underflow.

            // Get the correct schmidt values.
            // Remember that there are length+1 schmidt sites, since there is an extra "LC" somewhere on the chain.
            // We use the convention to take the schmidt values to the RIGHT of the current site.

            auto mps_self = state.get_mps_site(mps_site);
            long schmidt_host_site = mps_self.isCenter() ? mps_site : mps_site + 1;
            if(schmidt_host_site >=  state_len) prob = 0.0;
            else{
                tools::log->trace("Evaluating from {} ampl probability site {} | host {}", state.get_mps_site(mps_site).get_label(), mps_site, schmidt_host_site);
                auto &mps_host         = state.get_mps_site(schmidt_host_site);
                auto &schmidt_values   = mps_self.isCenter() ? mps_self.get_LC() : mps_host.get_L();
                eval_prob(schmidt_values);
            }
        }
    }

    void eval_from_B(const class_state_finite &state, long mps_site, const std::vector<Amplitude> &database) {
        // Start by calculating the mps sites that should be included in the amplitude
        // Remember that in eval_from_B we calculate the amplitude starting from the right-end of the chain
        if(not site or site.value() > mps_site) {
            // There are missing mps in the amplitude. Let's add them
            // First, check if it is already available in the amplitudes database
            // We can look it up by transforming the bits between 0 to bit_site to an integer
            // and use that as the index on amplitudes
            long state_len = state.get_length<long>();
            long state_pos = state.get_position<long>();
            long mps_rsite = state_len - 1 - mps_site;
            if(mps_site <= state_pos ) throw std::logic_error(fmt::format("eval_from_B: expected mps_site ({}) > state_pos ({})", mps_site,state_pos));

            if(mps_site < state_len-1) {
                // There is a chance to continue building an existing amplitude
                std::bitset<64> bits_index = 0;
                for(size_t i = 0; i < static_cast<size_t>(mps_rsite); i++) bits_index[i] = bits[i];
                const auto &c = database.at(bits_index.to_ulong()); // A cache item
                if(c.site and c.site.value() > mps_site) {
                    // Cache hit! No need to compute the amplitude from scratch
                    tools::log->trace("Cached   site {} | bits [{}] at idx {}", c.site.value(), c.to_string(), bits_index.to_ulong());
                    site = c.site;
                    ampl = c.ampl;
                }
            }

            if(ampl.size() == 0) { // Initialize
                ampl.resize(1);
                ampl.setConstant(1.0);
            }
            // Contract the missing mps
            for(const auto &mps : iter::reverse(state.mps_sites)) {
                long mps_pos = mps->get_position<long>();
                if(site and site.value() <= mps_pos) continue; // Fast-forward to the missing sites
                if(mps_site > mps_pos) break;                  // Contract up to the mps at mps_site
                if(ampl.size() != mps->get_chiR())
                    throw std::runtime_error(fmt::format("eval() failed for site {}: "
                                                         "mismatch in ampl({}) with size = {} and mps({}) with chiR = {} | bits {}",
                                                         mps_site, site.value(), ampl.size(), mps_pos, mps->get_chiL(), bits.to_string()));
                long                mps_rpos = state_len - 1 - mps->get_position<long>();
                long                size     = mps->get_chiL();
                std::array<long, 3> off      = {bits[static_cast<size_t>(mps_rpos)], 0, 0};
                std::array<long, 3> ext      = {1, mps->get_chiL(), mps->get_chiR()};

                // ampl never has a trailing Lambda, which means that we must SVD lambda out of B sites
                Eigen::Tensor<class_state_finite::Scalar, 1> temp =
                    mps->get_M_bare().slice(off, ext).contract(ampl, Textra::idx({2}, {0})).reshape(std::array<long, 1>{size});
                ampl = temp;
                // Update the current site
                site = mps_pos;
            }
            if(not site) throw std::logic_error(fmt::format("ampl has undefined site: should be {}", mps_site));
            if(site.value() != mps_site) throw std::logic_error(fmt::format("site ({}) != mps_site ({})",site.value(), mps_site ));
            // Get the correct schmidt values.
            // Remember that there are length+1 schmidt sites, since there is an extra "LC" somewhere on the chain.
            // We use the convention to take the schmidt values to the RIGHT of the current site.
            long  schmidt_host_site = mps_site - 1;
            if(schmidt_host_site < 0) prob = 0.0;
            else{
                auto &mps_host          = state.get_mps_site(schmidt_host_site);
                auto &schmidt_values = mps_host.isCenter() ? mps_host.get_LC() : mps_host.get_L();
                tools::log->trace("Evaluating from {} ampl probability site {} | host {}", state.get_mps_site(mps_site).get_label(), mps_site, schmidt_host_site);
                eval_prob(schmidt_values);
            }
        }
    }

    void eval(const class_state_finite &state, long mps_site, const std::vector<Amplitude> &database) {
        auto pos = state.get_position<long>();
        if(mps_site <= pos)
            eval_from_A(state, mps_site, database);
        else
            eval_from_B(state, mps_site, database);
    }
};

std::vector<double> compute_probability(const class_state_finite &state, long mps_pos, std::vector<Amplitude> &amplitudes, const std::vector<Amplitude> &database){
    auto state_pos = state.get_position<long>();
    auto state_len = state.get_length<long>();
    auto mps_rpos  = state_len - 1 - mps_pos;
    auto prob_size = mps_pos > state_pos ? mps_rpos + 2 : mps_pos+2;
    std::vector<std::optional<double>> probability(static_cast<size_t>(prob_size), std::nullopt);

    double p_sum = 0;
    for(auto &a : amplitudes) {
        if(probability[a.bits.count()]) {
//            if(tools::log->level() == spdlog::level::trace){
//                auto b = a; // Make a local copy
//                b.eval(state, mps_pos, database);
//                tools::log->info("Skipping P({:2}) = {:22.20f} | would have given {:22.20f}", b.bits.count(), probability[b.bits.count()].value(),
//                                 (b.prob ? b.prob.value() : std::numeric_limits<double>::quiet_NaN()));
//                if(b.prob and std::abs(probability[b.bits.count()].value() - b.prob.value()) > 1e-12)
//                    throw std::runtime_error(fmt::format("Skip would have missed significant contribution"));
//            }
            continue;
        }
        a.eval(state, mps_pos, database);

        if(a.prob) {
            std::string msg             = fmt::format("prob {:22.20f} to replace P({:2}) = {:22.20f}", a.prob.value(), a.bits.count(),
                                                      probability[a.bits.count()] ? probability[a.bits.count()].value() : std::numeric_limits<double>::quiet_NaN());
            probability[a.bits.count()] = a.prob;
            msg += fmt::format(" | p_sum {:.22f}", p_sum);
            p_sum = std::accumulate(probability.begin(), probability.end(), 0.0, [](auto &sum, auto &p) { return p ? sum + p.value() : sum; });
            msg += fmt::format(" -> {:.22f}", p_sum);
            tools::log->trace("{}", msg);
            if(std::abs(p_sum - 1.0) < 1e-16) break; // All particles found
        }
    }
    std::vector<double> p_vec;
    p_vec.reserve(probability.size());
    for(auto &p : probability) p_vec.emplace_back((p ? p.value() : 0.0));
    return p_vec;
}

std::vector<double> tools::finite::measure::number_entropies(class_state_finite &state) {
    auto save_pos = state.get_position<long>();
    auto save_dir = state.get_direction();
    tools::finite::mps::move_center_point_to_middle(state,state.find_largest_chi());
    tools::log->info("Labels: {}",state.get_labels());

    tools::log->enable_backtrace(200);
    auto state_pos = state.get_position<long>();
    auto state_len = state.get_length<long>();

    std::vector<double>    number_entropies(state.get_length<size_t>() + 1, 0.0); // Collects the resulting entropies
    std::vector<Amplitude> database;
    auto von_neumann_sum = [](double sum, const double p) {return p > 0 ? sum + p * std::log(p) : sum;};

    for(const auto &mps : state.mps_sites) {
        auto mps_pos = mps->get_position<long>();
        if(mps_pos > state_pos) break;
        if(mps->get_label() == "B") throw std::logic_error("Expected A/AC site, got B");
        // Generate a list of bit sequences of size prod_{i=0}^l spin_dim_i.
        // For spin-half this is just 2^l elements, where l is the mps position
        size_t num_bitseqs = std::accumulate(state.mps_sites.begin(), state.mps_sites.begin() + mps_pos + 1, 1ul, [](size_t &acc, const auto &mps) {
            if(mps->spin_dim() != 2) throw std::runtime_error("number_entropies: spin_dim() != 2 is not supported");
            return acc * mps->spin_dim();
        });

        std::vector<Amplitude> amplitudes;
        for(size_t count = 0; count < num_bitseqs; count++) amplitudes.emplace_back(Amplitude{count, std::nullopt, {}});

        auto probability = compute_probability(state,mps_pos,amplitudes,database);
        auto p_sum = std::accumulate(probability.begin(),probability.end(),0.0);
        if(std::abs(p_sum - 1.0) > 1e-8){
            tools::log->dump_backtrace();
            tools::log->info("P(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
            throw std::runtime_error(fmt::format("p_sum - 1.0 = {:.8e}", p_sum - 1.0));
        }
        tools::log->info("P(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
        double number_entropy = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
        number_entropies[mps->get_position<size_t>()] = std::abs(number_entropy);

        // Move the results into the database to use it in the next step
        database = amplitudes;
    }

    database.clear();

    for(const auto &mps : iter::reverse(state.mps_sites)) {
        auto mps_pos  = mps->get_position<long>();
        auto mps_rpos = state_len - 1 - mps_pos;
        if(mps_pos <= state_pos) break;
        if(mps->get_label() != "B") throw std::logic_error(fmt::format("Expected B site, got {}", mps->get_label()));

        // Generate a list of bit sequences of size prod_{i=0}^l spin_dim_i.
        // For spin-half this is just 2^l elements, where l is the mps position
        size_t num_bitseqs = std::accumulate(state.mps_sites.rbegin(), state.mps_sites.rbegin() + mps_rpos + 1, 1ul, [](size_t &acc, const auto &mps) {
            if(mps->spin_dim() != 2) throw std::runtime_error("number_entropies: spin_dim() != 2 is not supported");
            return acc * mps->spin_dim();
        });

        std::vector<Amplitude> amplitudes;
        for(size_t count = 0; count < num_bitseqs; count++) amplitudes.emplace_back(Amplitude{count, std::nullopt, {}});
        auto probability = compute_probability(state,mps_pos,amplitudes,database);
        auto p_sum = std::accumulate(probability.begin(),probability.end(),0.0);
        if(std::abs(p_sum - 1.0) > 1e-8){
            tools::log->dump_backtrace();
            tools::log->info("P(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
            throw std::runtime_error(fmt::format("p_sum - 1.0 = {:.8e}", p_sum - 1.0));
        }
        tools::log->info("P(n) = {:22.20f} = {:22.20f}", fmt::join(probability, ", "), p_sum);
        double number_entropy = -std::accumulate(probability.begin(), probability.end(), 0.0, von_neumann_sum);
        number_entropies[mps->get_position<size_t>()] = std::abs(number_entropy);

        // Move the results into the database to use it in the next step
        database = amplitudes;
    }
    tools::log->disable_backtrace();

    tools::finite::mps::move_center_point_to_pos_dir(state, save_pos,save_dir, state.find_largest_chi());
    return number_entropies;
}

