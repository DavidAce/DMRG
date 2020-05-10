//
// Created by david on 2019-01-29.
//

#include <general/nmspc_quantum_mechanics.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/debug.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>

void tools::finite::mps::initialize(class_state_finite &state, const std::string &model_type, size_t num_sites, size_t position) {
    log->info("Initializing mps with {} sites at position {}", num_sites, position);
    if(num_sites < 2) throw std::logic_error("Tried to initialize MPS with less than 2 sites");
    if(num_sites > 2048) throw std::logic_error("Tried to initialize MPS with more than 2048 sites");
    if(position >= num_sites) throw std::logic_error("Tried to initialize MPS at a position larger than the number of sites");

    size_t spin_dim = 2; // Default is a two-level system
    if(model_type == "ising_tf_rf_nn") spin_dim = settings::model::ising_tf_rf_nn::spin_dim;
    if(model_type == "ising_selfdual_tf_rf_nn") spin_dim = settings::model::ising_selfdual_tf_rf_nn::spin_dim;

    state.MPS_L.clear();
    state.MPS_R.clear();

    // Generate a simple MPS with all spins equal
    Eigen::Tensor<Scalar, 3> M(static_cast<long>(spin_dim), 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);
    M(0, 0, 0) = 0;
    M(1, 0, 0) = 1;
    L(0)       = 1;
    for(size_t site = 0; site < num_sites; site++) {
        if(site <= position)
            state.MPS_L.emplace_back(class_mps_site(M, L, site));
        else
            state.MPS_R.emplace_back(class_mps_site(M, L, site));
    }
    state.MPS_L.back().set_LC(L);
    if(state.MPS_L.size() + state.MPS_R.size() != num_sites) throw std::logic_error("Initialized MPS with wrong size");
    if(not state.get_MPS(position).isCenter()) throw std::logic_error("Initialized center matrix at the wrong position");
    if(state.get_position() != position) throw std::logic_error("Initialized MPS at the wrong position");
    state.site_update_tags = std::vector<bool>(num_sites, false);
}

void tools::finite::mps::random_product_state(class_state_finite &state, const std::string &parity_sector, const long state_number,
                                              const bool use_pauli_eigenstates)
/*!
 * There are many ways to random_product_state an initial product state state, based on the
 * arguments (parity_sector,state_number,use_pauli_eigenstates) = (string,long,true/false).
 * Let "+-sector" mean one of {"x","+x","-x","y","+y","-y", "z","+z","-z"}.

        a) ("+-sector"  ,+- ,t,f)   Set spinors to a random sequence of eigenvectors (up/down) of either
                                    sx, sy or sz pauli matrices (same pauli for all sites). If the global
                                    sign (+-) is omitted, a random sign is chosen with equal probabilities.
                                    In the x and z cases the full state will turn out to be entirely real,
                                    which improves performance.

        b) ("random"    ,+- ,f,f)   Set each spinor randomly on C2


        c) ("+-sector"  ,+- ,f,f)   Set each spinor randomly on C2 (i.e. case b) and then project the  full state
                                    to the given parity sector. If the global sign (+-) is omitted,  a random
                                    sign is chosen with equal probabilities. As a consequence of this, the
                                    full state will have always have nonzero imaginary part.

        d) ("randomAxis",+- ,f,f)   Randomly select one of {"x","y","z"} and go to case a).
        e) ("none"      ,+- ,f,f)   Does not random_product_state
        f) ("+-sector"  ,>=0,?,t)   Interpret seed_state as bitfield "01100010110..." and interpret these as
                                    up(0)/down(1) of either sx, sy or sz pauli matrices (same pauli for all sites)
 * Note: seed_state is only used if >= 0.
 * Note: we "use" the seed_state only once. Subsequent calls do not keep resetting the seed.
*/
{
    tools::log->debug("Randomizing mps into sector {}", parity_sector);
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_have_been_updated(false);

    if(state_number >= 0)
        internals::set_product_state_in_parity_sector_from_bitset(state, parity_sector, state_number);
    else
        internals::set_product_state_randomly(state, parity_sector, use_pauli_eigenstates);
    tools::finite::mps::rebuild_environments(state);
}

void tools::finite::mps::random_current_state(class_state_finite &state, const std::string &parity_sector1, const std::string &parity_sector2) {
    Eigen::MatrixXcd paulimatrix1;
    Eigen::MatrixXcd paulimatrix2;
    if(parity_sector1 == "x")
        paulimatrix1 = qm::spinOneHalf::sx;
    else if(parity_sector1 == "y")
        paulimatrix1 = qm::spinOneHalf::sy;
    else if(parity_sector1 == "z")
        paulimatrix1 = qm::spinOneHalf::sz;
    else
        paulimatrix1 = qm::spinOneHalf::Id;
    if(parity_sector2 == "x")
        paulimatrix2 = qm::spinOneHalf::sx;
    else if(parity_sector2 == "y")
        paulimatrix2 = qm::spinOneHalf::sy;
    else if(parity_sector2 == "z")
        paulimatrix2 = qm::spinOneHalf::sz;
    else
        paulimatrix2 = qm::spinOneHalf::Id;
    //    auto [mpos,L,R] = qm::mpo::random_pauli_mpos(paulimatrix,state.get_length());
    auto chi_lim      = state.find_largest_chi();
    auto [mpos, L, R] = qm::mpo::random_pauli_mpos_x2(paulimatrix1, paulimatrix2, state.get_length());
    tools::finite::ops::apply_mpos(state, mpos, L, R);
    tools::finite::mps::normalize(state, chi_lim);
    tools::finite::debug::check_integrity(state);
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, "x");
}

void tools::finite::mps::rebuild_environments(class_state_finite &state) {
    tools::log->trace("Rebuilding environments");
    if(state.MPS_L.size() != state.MPO_L.size())
        throw std::runtime_error(fmt::format("Size mismatch in MPSL and MPOL: {} != {}", state.MPS_L.size(), state.MPO_L.size()));
    if(state.MPS_R.size() != state.MPO_R.size())
        throw std::runtime_error(fmt::format("Size mismatch in MPSR and MPOR: {} != {}", state.MPS_R.size(), state.MPO_R.size()));
    // Generate new environments

    {
        state.ENV_L.clear();
        state.ENV2_L.clear();

        auto mpsL_it = state.MPS_L.begin();
        auto mpoL_it = state.MPO_L.begin();
        auto ENV_L   = class_environment("L", *mpsL_it, *mpoL_it->get());     // Initialized envs
        auto ENV2_L  = class_environment_var("L", *mpsL_it, *mpoL_it->get()); // Initialized envs
        while(mpsL_it != state.MPS_L.end() and mpoL_it != state.MPO_L.end()) {
            if(ENV_L.hasNaN()) throw std::runtime_error("ENV_L " + std::to_string(ENV_L.get_position()) + " has NAN's");
            state.ENV_L.emplace_back(ENV_L);
            state.ENV2_L.emplace_back(ENV2_L);
            if(mpsL_it->get_position() != state.ENV_L.back().get_position())
                throw std::runtime_error(fmt::format("Size mismatch in MPSL and ENVL: {} != {}", mpsL_it->get_position(), state.ENV_L.back().get_position()));
            if(mpsL_it->get_chiL() != state.ENV_L.back().block.dimension(0))
                throw std::runtime_error(
                    fmt::format("Size mismatch in MPSL and ENVL dimensions {} != {}", mpsL_it->get_chiL(), state.ENV_L.back().block.dimension(2)));

            ENV_L  = ENV_L.enlarge(*mpsL_it, *mpoL_it->get());
            ENV2_L = ENV2_L.enlarge(*mpsL_it, *mpoL_it->get());
            mpsL_it++;
            mpoL_it++;
        }
    }

    {
        state.ENV_R.clear();
        state.ENV2_R.clear();

        auto mpsR_it = state.MPS_R.rbegin();
        auto mpoR_it = state.MPO_R.rbegin();
        auto ENV_R   = class_environment("R", *mpsR_it, *mpoR_it->get());
        auto ENV2_R  = class_environment_var("R", *mpsR_it, *mpoR_it->get());
        while(mpsR_it != state.MPS_R.rend() and mpoR_it != state.MPO_R.rend()) {
            if(ENV_R.hasNaN()) throw std::runtime_error("ENV_R " + std::to_string(ENV_R.get_position()) + " has NAN's");
            state.ENV_R.emplace_front(ENV_R);
            state.ENV2_R.emplace_front(ENV2_R);
            if(mpsR_it->get_position() != state.ENV_R.front().get_position())
                throw std::runtime_error(fmt::format("Size mismatch in MPSR and ENVR: {} != {}", mpsR_it->get_position(), state.ENV_R.front().get_position()));
            if(mpsR_it->get_chiR() != state.ENV_R.front().block.dimension(0))
                throw std::runtime_error(
                    fmt::format("Size mismatch in MPSR and ENVR dimensions {} != {}", mpsR_it->get_chiR(), state.ENV_R.front().block.dimension(2)));
            ENV_R  = ENV_R.enlarge(*mpsR_it, *mpoR_it->get());
            ENV2_R = ENV2_R.enlarge(*mpsR_it, *mpoR_it->get());
            mpsR_it++;
            mpoR_it++;
        }
    }
}

void tools::finite::mps::project_to_closest_parity_sector(class_state_finite &state, std::string parity_sector) {
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, parity_sector);
}
