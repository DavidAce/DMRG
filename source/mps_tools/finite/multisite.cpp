//
// Created by david on 2019-06-24.
//


#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <sim_parameters/nmspc_sim_settings.h>

std::vector<size_t> mpstools::finite::multisite::compute_best_jump(class_finite_chain_state &state,  mpstools::finite::opt::OptSpace optSpace){
    using namespace Textra;
    size_t position  = state.get_position();
    size_t direction = state.get_direction();
    size_t length    = state.get_length();
    std::vector<long>   costs;
    std::vector<size_t> sites;
    std::vector<Eigen::DSizes<long,3>> dims;

    while(position >= 0 and position < length){
        dims.emplace_back(state.get_G(position).dimensions());
        long cost = direction > 1 ? dims.back()[1]*dims.front()[2] : dims.front()[1]*dims.back()[2]  ;
        for (auto &d : dims ){cost *= d[0];}
        costs.push_back(cost);
        sites.push_back(position);
        position += direction;
    }
    std::reverse(costs.begin(),costs.end());

    std::cout << "JUMP COSTS \n" << costs << std::endl;

    // Evaluate best cost. Threshold depends on optSpace
    // Case 1: All costs are equal              -> take all sites
    // Case 2: Costs increase indefinitely      -> take until threshold
    // Case 3: Costs increase and saturate      -> take until threshold
    long threshold = 0;
    switch(optSpace){
        case mpstools::finite::opt::OptSpace::FULL    : threshold = 2 * 2 * 16 * 16; break;
        case mpstools::finite::opt::OptSpace::PARTIAL : threshold = 2 * 2 * 32 * 16; break;
        case mpstools::finite::opt::OptSpace::DIRECT  : threshold = 2 * 2 * settings::xdmrg::chi_max *  settings::xdmrg::chi_max; break;
    }

    auto costsmap = Eigen::Map<Eigen::Array<long, Eigen::Dynamic,1>>(costs.data(),costs.size());
    bool allequal = (costsmap == costsmap(0)).all();
    for (auto & c : costs){
        if (allequal){std::cout << "allequal\n"; break;}
        if (sites.size() == 1){break;}
        if (sites.size() == 0){throw std::logic_error("No sites for a jump");}
        if (c <= threshold){break;}
        else{
            sites.pop_back();
        }
    }
    if (direction == -1){std::reverse(sites.begin(),sites.end());}
    std::cout << "SITES \n" << sites << std::endl;
    return sites;

}
