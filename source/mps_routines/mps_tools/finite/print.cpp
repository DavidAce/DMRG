//
// Created by david on 2019-01-29.
//


#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_environment.h>
#include <general/nmspc_tensor_extra.h>
#include <string>
#include <iomanip>


using Scalar         = std::complex<double>;

void MPS_Tools::Finite::Print::print_full_state(const class_finite_chain_state &state) {
    for (auto & mps : state.get_MPS_L()){
        std::cout << "MPS " << mps.get_position() << "  :\n";
        std::cout << "  L:\n"<< mps.get_L() << '\n';
        std::cout << "  G:\n"<< mps.get_G() << '\n';
    }
    std::cout << "  LC:\n"<< state.get_MPS_C();
    for (auto & mps : state.get_MPS_R()){
        std::cout << "MPS " << mps.get_position() << "  :\n";
        std::cout << "  G:\n"<< mps.get_G() << '\n';
        std::cout << "  L:\n"<< mps.get_L() << '\n';
    }
}



void MPS_Tools::Finite::Print::print_state(const class_finite_chain_state &state){
    auto & MPS_L  = state.get_MPS_L();
    auto & MPS_R  = state.get_MPS_R();
    auto & MPS_C  = state.get_MPS_C();
    auto & ENV_L  = state.get_ENV_L();
    auto & ENV_R  = state.get_ENV_R();

    int i = 0;
    std::cout << std::setprecision(10);
    std::cout << "State length              : " << state.get_length() << std::endl;
    std::cout << "State position            : "    << state.get_position() << std::endl;
    std::cout << "Environment L size        : "    << ENV_L.size() << std::endl;
    std::cout << "Environment R size        : "    << ENV_R.size() << std::endl;

    if(!ENV_L.empty()){std::cout << "ENV_L[" << std::setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().size << "  <--- Also current environment L" << std::endl;}

    for(auto &it : MPS_L){
        std::cout << "L[" << std::setw(3) << i  <<  "]: " << it.get_L().dimensions()  << " "
                  << "G[" << std::setw(3) << i  <<  "]: " << it.get_G().dimensions()  << " ";
        if(&it == &MPS_L.back()){
            std::cout << " <--- Position A";
        }
        std::cout << std::endl;
        i++;
    }
    std::cout << "L[" << std::setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    for(auto &it : MPS_R){
        std::cout << "G[" << std::setw(3) << i  <<  "]: " << it.get_G().dimensions()  << " "
                  << "L[" << std::setw(3) << i  <<  "]: " << it.get_L().dimensions()  << " ";
        if(&it == &MPS_R.front()){
            std::cout << " <--- Position B" ;
        }
        std::cout << std::endl;
        i++;
    }
    if(!ENV_R.empty()){std::cout << "ENV_R[" << std::setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}

}


void MPS_Tools::Finite::Print::print_state_compact(const class_finite_chain_state &state){
    auto & MPS_L  = state.get_MPS_L();
    auto & MPS_R  = state.get_MPS_R();
    auto & MPS_C  = state.get_MPS_C();
    auto & ENV_L  = state.get_ENV_L();
    auto & ENV_R  = state.get_ENV_R();

    std::cout << std::setprecision(10);

    std::cout << "State length              : "    << state.get_length() << std::endl;
    std::cout << "State position            : "    << state.get_position() << std::endl;
    std::cout << "Environment L size        : "    << ENV_L.size() << std::endl;
    std::cout << "Environment R size        : "    << ENV_R.size() << std::endl;
    if(!ENV_L.empty()){std::cout << "ENV_L[" <<std::setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().size << "  <--- Also current environment L" << std::endl;}
    if(!MPS_L.empty()){std::cout << "MPS_L[" <<std::setw(3) << MPS_L.size()-1 << "]: " << MPS_L.back().get_G().dimensions() <<  "   <--- Also current iteration A" << std::endl;}
    std::cout << "L[" << std::setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    if(!MPS_R.empty()){std::cout << "MPS_R[" <<std::setw(3) << MPS_R.size()-1 << "]: " << MPS_R.front().get_G().dimensions() << "   <--- Also current iteration B" << std::endl;}
    if(!ENV_R.empty()){std::cout << "ENV_R[" <<std::setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}
}




void MPS_Tools::Finite::Print::print_hamiltonians(const class_finite_chain_state &state) {
    auto & MPO_L  = state.get_MPO_L();
    auto & MPO_R  = state.get_MPO_R();

    MPO_L.begin()->get()->print_parameter_names();
    for(auto &it : MPO_L){
        it->print_parameter_values();
    }
    for(auto &it : MPO_R){
        it->print_parameter_values();
    }
}
