//
// Created by david on 2019-01-29.
//

#include <iomanip>
#include <state/class_state_finite.h>
#include <string>
#include <tools/common/log.h>
#include <tools/finite/print.h>

using Scalar         = std::complex<double>;

void tools::finite::print::print_full_state(const class_state_finite &state) {
    
    for (auto & mps : state.MPS_L){
        std::cout << "MPS " << mps.get_position() << "  :\n";
        std::cout << "  L:\n"<< mps.get_L() << '\n';
        std::cout << "  M:\n"<< mps.get_M() << '\n';
    }
    std::cout << "  LC:\n"<< state.MPS_L.back().get_LC() << '\n';
    for (auto & mps : state.MPS_R){
        std::cout << "MPS " << mps.get_position() << "  :\n";
        std::cout << "  M:\n"<< mps.get_M() << '\n';
        std::cout << "  L:\n"<< mps.get_L() << '\n';
    }
}



void tools::finite::print::print_state(const class_state_finite &state){
    using namespace Textra;
    tools::log->info("State length              {}", state.get_length());
    tools::log->info("State position            {}", state.get_position());
    tools::log->info("Environment L size        {}", state.ENV_L.size());
    tools::log->info("Environment R size        {}", state.ENV_R.size());

    for(size_t pos = 0; pos < state.get_length(); pos++){
        std::string tag;
        if(pos == state.get_position())   tag = "<---- Position A";
        if(pos == state.get_position()+1) tag = "<---- Position B";
        const auto & mps = state.get_MPS(pos);
        if(pos <= state.get_position()) {
            const auto & env = state.get_ENVL(pos).block;
            tools::log->info("Pos {:2}: L [{:^4}] M [{:>2} {:>3} {:>3}] ENV [{:>3} {:>3} {:>2}] {}", pos, mps.get_L().dimension(0), mps.get_spin_dim(),
                             mps.get_chiL(), mps.get_chiR(), env.dimension(0), env.dimension(1), env.dimension(2), tag);
        } else {
            const auto & env = state.get_ENVR(pos).block;
            tools::log->info("Pos {:2}: M [{:>2} {:>3} {:>3}] L [{:^4}] ENV [{:>3} {:>3} {:>2}] {}", pos, mps.get_spin_dim(), mps.get_chiL(), mps.get_chiR(),
                             mps.get_L().dimension(0), env.dimension(0), env.dimension(1), env.dimension(2), tag);
        }
        if(state.get_MPS(pos).isCenter())
        tools::log->info("Pos {:2}: L [{:^4}] {:>44}", pos, state.get_MPS(pos).get_L().dimension(0),"<---- Center");
    }
}


void tools::finite::print::print_state_compact(const class_state_finite &state){
    using namespace Textra;
    auto & MPS_L  = state.MPS_L;
    auto & MPS_R  = state.MPS_R;
    auto & MPS_C  = state.midchain_bond();
    auto & ENV_L  = state.ENV_L;
    auto & ENV_R  = state.ENV_R;

    std::cout << std::setprecision(10);

    std::cout << "State length              : "    << state.get_length() << std::endl;
    std::cout << "State position            : "    << state.get_position() << std::endl;
    std::cout << "Environment L size        : "    << ENV_L.size() << std::endl;
    std::cout << "Environment R size        : "    << ENV_R.size() << std::endl;
    if(!ENV_L.empty()){std::cout << "ENV_L[" <<std::setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().sites << "  <--- Also current environment L" << std::endl;}
    if(!MPS_L.empty()){std::cout << "MPS_L[" <<std::setw(3) << MPS_L.size()-1 << "]: " << MPS_L.back().get_M().dimensions() <<  "   <--- Also current M" << std::endl;}
    std::cout << "L[" << std::setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    if(!MPS_R.empty()){std::cout << "MPS_R[" <<std::setw(3) << MPS_R.size()-1 << "]: " << MPS_R.front().get_M().dimensions() << "   <--- Also current M" << std::endl;}
    if(!ENV_R.empty()){std::cout << "ENV_R[" <<std::setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().sites << " <--- Also current environment R"  << std::endl;}
}




void tools::finite::print::print_hamiltonians(const class_state_finite &state) {
    auto & MPO_L  = state.MPO_L;
    auto & MPO_R  = state.MPO_R;
    if (MPO_L.empty()) throw std::runtime_error("MPO_L is empty. Can't print hamiltonian");
    if (MPO_R.empty()) throw std::runtime_error("MPO_R is empty. Can't print hamiltonian");

    MPO_L.begin()->get()->print_parameter_names();
    for(auto &it : MPO_L){
        it->print_parameter_values();
    }
    for(auto &it : MPO_R){
        it->print_parameter_values();
    }
}
