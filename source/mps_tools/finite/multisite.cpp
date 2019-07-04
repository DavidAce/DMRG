//
// Created by david on 2019-06-24.
//


#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <sim_parameters/nmspc_sim_settings.h>

std::list<size_t> mpstools::finite::multisite::generate_site_list(class_finite_chain_state &state, long threshold){
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
        long spindim = 1;
        for (auto &d : dims ){spindim *= d[0];}
        cost *= spindim;
        std::cout << "position: " << position << " total dims: [ " << spindim << " " << chiL << " " << chiR << " ] dims: " << dims.back() << " cost: " << cost << std::endl;
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

    auto costsmap = Eigen::Map<Eigen::Array<long, Eigen::Dynamic,1>>(costs.data(),costs.size());
    for (auto & c : costs){
//        if (allequal){std::cout << "allequal\n"; break;}
        if (sites.size() <= 2){std::cout << "at least two sites kept \n"; break;}
        else if (sites.empty()){throw std::logic_error("No sites for a jump");}
        else if (c <= threshold and sites.size() <= 8){std::cout << "good threshold found: " << c << '\n';break;}
        else{
            sites.pop_back();
        }
    }
    if (direction == -1){std::reverse(sites.begin(),sites.end());}
    std::cout << "SITES \n" << sites << std::endl;
    return sites;

}
