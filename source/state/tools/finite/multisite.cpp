//
// Created by david on 2019-06-24.
//


#include <state/tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <simulation/nmspc_settings.h>
#include <spdlog/fmt/bundled/ranges.h>

std::list<size_t> tools::finite::multisite::generate_site_list(class_finite_state &state, long threshold){
    using namespace Textra;
    int    direction = state.get_direction();
    size_t position  = state.get_position();
    size_t length    = state.get_length();
    std::vector<long> costs;
    std::list<size_t> sites;
    std::vector<Eigen::DSizes<long,3>> dims;
    if (direction == -1)position++; // If going to the right, take position to be the site on the right of the center bond.
    while(position >= 0 and position < length){
        dims.emplace_back(state.get_G(position).dimensions());
        long chiL = direction == 1 ? dims.front()[1] : dims.back() [1];
        long chiR = direction == 1 ? dims.back() [2] : dims.front()[2];
        long cost = chiL * chiR;
        for (auto &d : dims ){cost *= d[0];}
        costs.push_back(cost);
        sites.push_back(position);
        position += direction;
    }
    std::reverse(costs.begin(),costs.end());


    // Evaluate best cost. Threshold depends on optSpace
    // Case 1: All costs are equal              -> take all sites
    // Case 2: Costs increase indefinitely      -> take until threshold
    // Case 3: Costs increase and saturate      -> take until threshold
    auto costsmap = Eigen::Map<Eigen::Array<long, Eigen::Dynamic,1>>(costs.data(),costs.size());
    bool allequal = (costsmap == costsmap(0)).all();
    std::string reason;
    for (auto & c : costs){
        using namespace settings::precision;
        if (allequal)                                   {reason = "equal costs: " + std::to_string(c); break;}
        else if (sites.size() <= 2)                     {reason = "at least two sites were kept"; break;}
        else if (c <= threshold
                and sites.size() <= MaxSitesMultiDmrg)  {reason = "good threshold found: " + std::to_string(c) ;break;}
        else if (sites.empty())                         {throw std::logic_error("No sites for a jump");}
        else{
            sites.pop_back();
        }
    }
    tools::log->debug("Problem sizes: {}", costs);
    tools::log->debug("Chosen sites {}. Reason: {}", sites, reason);
    if (direction == -1){std::reverse(sites.begin(),sites.end());}
    state.active_sites = sites;
    return sites;
}


using namespace Textra;
using Scalar = class_finite_state::Scalar;


double tools::finite::measure::multisite::energy(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
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
        throw std::runtime_error("Energy has an imaginary part: " + std::to_string(std::real(E(0))) + " + i " + std::to_string(std::imag(E(0))));
    }
    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");
    return std::real(E(0));
}


double tools::finite::measure::multisite::energy_per_site(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
        return multisite::energy(state,multitheta)/state.get_length();
}


double tools::finite::measure::multisite::energy_variance(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
    double energy = tools::finite::measure::multisite::energy(state, multitheta);
    auto multimpo   = state.get_multimpo();
    auto & env2L    = state.get_ENV2L(state.active_sites.front()).block;
    auto & env2R    = state.get_ENV2R(state.active_sites.back()).block;
    Eigen::Tensor<Scalar, 0> H2 =
            env2L
            .contract(multitheta                 , idx({0}  ,{1}))
            .contract(multimpo                   , idx({3,1},{2,0}))
            .contract(multimpo                   , idx({4,1},{2,0}))
            .contract(multitheta.conjugate()     , idx({4,0},{0,1}))
            .contract(env2R                      , idx({0,3,1,2},{0,1,2,3}));
    return std::abs(H2(0) - energy*energy);
}


double tools::finite::measure::multisite::energy_variance_per_site(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
        return multisite::energy_variance(state,multitheta)/state.get_length();
}


