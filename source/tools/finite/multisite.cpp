//
// Created by david on 2019-06-24.
//



#ifdef _OPENMP
#include <omp.h>
#define EIGEN_USE_THREADS
#include <unsupported/Eigen/CXX11/Tensor>
Eigen::ThreadPool       tp_multisite (Eigen::nbThreads());
Eigen::ThreadPoolDevice dev_multisite(&tp,Eigen::nbThreads());
#else
#include <unsupported/Eigen/CXX11/Tensor>
Eigen::DefaultDevice dev_multisite;
#endif

#include <tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <simulation/nmspc_settings.h>
#include <spdlog/fmt/bundled/ranges.h>


Eigen::DSizes<long,3> tools::finite::multisite::get_dimensions  (const class_finite_state &state, const std::list<size_t> &list_of_sites){
    if (list_of_sites.empty()) return  Eigen::DSizes<long,3>{0,0,0};
    Eigen::DSizes<long,3> dimensions;
    int direction = list_of_sites.back() >= list_of_sites.front() ? 1 : -1;
    if (direction == 1){
        dimensions[1] = state.get_G(list_of_sites.front()).dimension(1);
        dimensions[2] = state.get_G(list_of_sites.back()) .dimension(2);
    }
    else{
        dimensions[1] = state.get_G(list_of_sites.back()) .dimension(1);
        dimensions[2] = state.get_G(list_of_sites.front()).dimension(2);
    }

    dimensions[0] = 1;
    for (auto & site : list_of_sites){
        dimensions[0] *= state.get_G(site).dimension(0);
    }
    return dimensions;
}


size_t tools::finite::multisite::get_problem_size(const class_finite_state &state, const std::list<size_t> &list_of_sites){
    auto dims = get_dimensions(state,list_of_sites);
    return dims[0]*dims[1]*dims[2];
}




std::list<size_t> tools::finite::multisite::generate_site_list(class_finite_state &state, const size_t threshold, const size_t max_sites){
    tools::log->trace("Activating sites. Threshold = {}, Max sites = {}", threshold,max_sites);
    using namespace Textra;
    int    direction = state.get_direction();
    size_t position  = state.get_position();
    size_t length    = state.get_length();
    if (direction == -1)position++; // If going to the left, take position to be the site on the right of the center bond.
    std::list<size_t> costs;
    std::list<size_t> sites;
    std::vector<Eigen::DSizes<long,3>> dims;
    while(position >= 0 and position < length){
        sites.emplace_back(position);
        costs.emplace_back(get_problem_size(state,sites));
        position += direction;
    }
    tools::log->debug("Activation problem sizes: {}", costs);
    // Evaluate best cost. Threshold depends on optSpace
    // Case 1: All costs are equal              -> take all sites
    // Case 2: Costs increase indefinitely      -> take until threshold
    // Case 3: Costs increase and saturate      -> take until threshold

    std::string reason;
    while (true){
        bool allequal = std::all_of(costs.begin(), costs.end(), [costs](size_t c) { return c == costs.front(); });
        auto c = costs.back();
        if (c <= threshold  and sites.size() <= max_sites) {reason = "good threshold found: " + std::to_string(c) ;break;}
        else if (sites.size() <= 2)                        {reason = "at least two sites were kept"; break;}
        else if (allequal and sites.size() <= max_sites)   {reason = "equal costs: " + std::to_string(c); break;}
        else if (sites.size() == 1)                        {throw std::logic_error("At least two sites required!");}
        else if (sites.empty())                            {throw std::logic_error("No sites for a jump");}
        else{
            sites.pop_back();
            costs.pop_back();
        }
    }
    if (direction == -1){std::reverse(sites.begin(),sites.end());}
    tools::log->debug("Chosen sites {}. Reason: {}", sites, reason);
    state.active_sites = sites;
    return sites;
}


using namespace Textra;
using Scalar = class_finite_state::Scalar;


double tools::finite::measure::multisite::internal::significant_digits(double H2, double E2){
    double max_digits    = std::numeric_limits<double>::max_digits10;
    double lost_bits     = -std::log2(1.0 - std::abs(std::min(H2,E2)/std::max(H2,E2)));
    double lost_digits   = std::log10(std::pow(2.0,lost_bits));
//    tools::log->trace("Significant digits: {}",std::floor(max_digits - lost_digits));
    return digits = std::floor(max_digits - lost_digits);
}

double tools::finite::measure::multisite::energy_minus_energy_reduced(const class_finite_state &state, const Eigen::Tensor<Scalar,3> & multitheta){
    tools::common::profile::t_ene.tic();
    auto multimpo   = state.get_multimpo();
    auto & envL     = state.get_ENVL(state.active_sites.front()).block;
    auto & envR     = state.get_ENVR(state.active_sites.back()).block;

    Eigen::Tensor<Scalar, 0>  E =
            envL
                    .contract(multitheta,                               idx({0},{1}))
                    .contract(multimpo,                                 idx({2,1},{2,0}))
                    .contract(multitheta.conjugate(),                   idx({3,0},{0,1}))
                    .contract(envR,                                     idx({0,2,1},{0,1,2}));
    if(abs(imag(E(0))) > 1e-10 ){
        tools::log->critical(fmt::format("Energy has an imaginary part: {:.16f} + i {:.16f}",std::real(E(0)), std::imag(E(0))));
//        throw std::runtime_error("Energy has an imaginary part: " + std::to_string(std::real(E(0))) + " + i " + std::to_string(std::imag(E(0))));
    }

    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");
    double ene = std::real(E(0));
    if (std::isnan(ene) or std::isinf(ene)) throw std::runtime_error(fmt::format("Energy is invalid: {}", ene));
    tools::common::profile::t_ene.toc();
    return  ene;
}


double tools::finite::measure::multisite::energy(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
    // We want to measure energy accurately always.
    // Since the state can be reduced, the true energy is always
    // E + E_reduced
//    double e_minus_e_reduced = multisite::energy_minus_energy_reduced(state,multitheta);
//    double e_reduced = state.get_energy_reduced();
//    tools::log->debug("Energy minus Energy_reduced = {}",e_minus_e_reduced);
//    tools::log->debug("Energy_reduced              = {}",e_reduced);
    return multisite::energy_minus_energy_reduced(state,multitheta) + state.get_energy_reduced();
}


double tools::finite::measure::multisite::energy_per_site(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
        return multisite::energy(state,multitheta)/state.get_length();
}


double tools::finite::measure::multisite::energy_variance(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
    // Depending on whether the state is reduced or not we get different formulas.
    // Luckily, the variance is independent of offsets.
    // If the state is not reduced we get Var H = H^2 - E^2 =  H2 - energy*energy
    // IF the state is reduced we get Var H = (H-E_red) - (E-E_red)^2 = H2 - energy_minus_energy_reduced^2

    tools::common::profile::t_var.tic();
    auto multimpo   = state.get_multimpo();
    auto & env2L    = state.get_ENV2L(state.active_sites.front()).block;
    auto & env2R    = state.get_ENV2R(state.active_sites.back()).block;


    auto dsizes      = state.active_dimensions();
    size_t log2chiL  = std::log2(dsizes[1]);
    size_t log2chiR  = std::log2(dsizes[2]);
    size_t log2spin  = std::log2(dsizes[0]);
    Eigen::Tensor<Scalar, 0> H2;
    if (log2spin > log2chiL + log2chiR){
        if (log2chiL > log2chiR){
//            tools::log->trace("H2 path: log2spin > log2chiL + log2chiR  and  log2chiL > log2chiR ");
            Eigen::Tensor<Scalar,3> theta = multitheta.shuffle(Textra::array3{1,0,2});
            H2.device(dev_multisite) =
                    theta
                            .contract(env2L              , Textra::idx({0}, {0}))
                            .contract(multimpo           , Textra::idx({0,3}, {2,0}))
                            .contract(env2R              , Textra::idx({0,3}, {0,2}))
                            .contract(multimpo           , Textra::idx({2,1,4}, {2,0,1}))
                            .contract(theta.conjugate()  , Textra::idx({2,0,1}, {1,0,2}));
        }

        else{
//            tools::log->trace("H2 path: log2spin > log2chiL + log2chiR  and  log2chiL <= log2chiR ");
            Eigen::Tensor<Scalar,3> theta = multitheta.shuffle(Textra::array3{2,0,1});
            H2.device(dev_multisite) =
                    theta
                            .contract(env2R              , Textra::idx({0}, {0}))
                            .contract(multimpo           , Textra::idx({0,3}, {2,1}))
                            .contract(env2L              , Textra::idx({0,3}, {0,2}))
                            .contract(multimpo           , Textra::idx({2,4,1}, {2,0,1}))
                            .contract(theta.conjugate()  , Textra::idx({2,1,0}, {1,2,0}));
        }

    }else{
//        tools::log->trace("H2 path: log2spin <= log2chiL + log2chiR");
        Eigen::Tensor<Scalar,3> theta = multitheta.shuffle(Textra::array3{1,0,2});
        H2.device(dev_multisite) =
                theta
                        .contract(env2L              , Textra::idx({0}, {0}))
                        .contract(multimpo           , Textra::idx({0,3}, {2,0}))
                        .contract(multimpo           , Textra::idx({4,2}, {2,0}))
                        .contract(env2R              , Textra::idx({0,2,3}, {0,2,3}))
                        .contract(theta.conjugate()  , Textra::idx({1,0,2}, {1,0,2}));
    }




//
//
//
//    Eigen::Tensor<Scalar, 0> H2 =
//            env2L
//            .contract(multitheta                 , idx({0}  ,{1}))
//            .contract(multimpo                   , idx({3,1},{2,0}))
//            .contract(multimpo                   , idx({4,1},{2,0}))
//            .contract(multitheta.conjugate()     , idx({4,0},{0,1}))
//            .contract(env2R                      , idx({0,3,1,2},{0,1,2,3}));
    tools::common::profile::t_var.toc();
    double energy;
    if (state.isReduced()){
        energy = multisite::energy_minus_energy_reduced(state,multitheta);
    }else{
        energy = multisite::energy(state, multitheta);
    }
    double E2 = energy * energy;
    double var = std::abs(H2(0) - E2);
    if (std::isnan(var) or std::isinf(var)) throw std::runtime_error(fmt::format("Variance is invalid: {}", var));
    internal::significant_digits(std::abs(H2(0)),E2);
    return var;
}


double tools::finite::measure::multisite::energy_variance_per_site(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
        return multisite::energy_variance(state,multitheta)/state.get_length();
}



double tools::finite::measure::multisite::energy(const class_finite_state &state){
    if (state.measurements.energy)  return state.measurements.energy.value();
    if (state.active_sites.empty()) return tools::finite::measure::energy(state);
    tools::common::profile::t_ene.tic();
    auto theta = state.get_multitheta();
    tools::common::profile::t_ene.toc();
    state.measurements.energy = multisite::energy(state,theta);
    return state.measurements.energy.value();
}

double tools::finite::measure::multisite::energy_per_site(const class_finite_state &state){
    if (state.measurements.energy_per_site){return state.measurements.energy_per_site.value();}
    else{
        if (state.active_sites.empty()) return tools::finite::measure::energy_per_site(state);
        state.measurements.energy_per_site = multisite::energy(state)/state.get_length();
        return state.measurements.energy_per_site.value();
    }
}

double tools::finite::measure::multisite::energy_variance(const class_finite_state &state){
    if (state.measurements.energy_variance_mpo){return state.measurements.energy_variance_mpo.value();}
    else{
        if (state.active_sites.empty()) return tools::finite::measure::energy_variance(state);
        tools::common::profile::t_var.tic();
        auto theta = state.get_multitheta();
        tools::common::profile::t_var.toc();
        state.measurements.energy_variance_mpo = multisite::energy_variance(state,theta);
        return state.measurements.energy_variance_mpo.value();
    }}

double tools::finite::measure::multisite::energy_variance_per_site(const class_finite_state &state){
    if (state.measurements.energy_variance_per_site){return state.measurements.energy_variance_per_site.value();}
    else{
        if (state.active_sites.empty()) return tools::finite::measure::energy_variance_per_site(state);
        state.measurements.energy_variance_per_site = multisite::energy_variance(state)/state.get_length();
        return state.measurements.energy_variance_per_site.value();
    }
}




