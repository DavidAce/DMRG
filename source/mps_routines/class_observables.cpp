//
// Created by david on 2017-11-12.
//

#include <iomanip>
#include "class_observables.h"
#include "general/class_eig_spectra_wrapper.h"
class_observables::class_observables(class_superblock &superblockRef, SimulationType sim_)
                    :superblock(superblockRef),
                    sim(sim_)
{
    Sim2String = {
            {SimulationType::iDMRG,     "iDMRG: "},
            {SimulationType::fDMRG,     "fDMRG: "},
            {SimulationType::FES_iDMRG, "FES_iDMRG: "},
            {SimulationType::iTEBD,     "iTEBD: "},
            {SimulationType::FES_iTEBD, "FES_iTEBD: "}
    };
}

double class_observables::get_expectationvalue(const Textra::Tensor<4,Scalar> &MPO){
    Textra::Tensor<4,Scalar> theta  = superblock.MPS.get_theta();
    Textra::Tensor<0,Scalar> result =  superblock.Lblock.block
            .contract(theta.conjugate(),idx<1>({0},{2}))
            .contract(MPO,              idx<2>({1,2},{0,2}))
            .contract(MPO,              idx<2>({3,1},{0,2}))
            .contract(theta,            idx<3>({0,2,4},{2,0,1}))
            .contract(superblock.Rblock.block,     idx<3>({0,1,2},{0,2,1}));
    return static_cast<double>(result(0))/ superblock.chain_length;
}

double class_observables::get_energy(){
    if (Sim2String[sim].find("TEBD") != std::string::npos ){
        Textra::Tensor<4,Scalar> theta  = superblock.MPS.get_theta();
        Textra::Tensor<0,Scalar> E1     =  superblock.H.asTensor
                                           .contract(theta.conjugate(), idx<2>({0,1},{0,1}))
                                           .contract(theta,             idx<4>({0,1,2,3},{0,1,2,3}));
        return static_cast<double>(E1(0));
    }else {
        return get_expectationvalue(superblock.H.M);
    }
}

double class_observables::get_entropy(){
    Textra::Tensor<0,Scalar> result1  = -superblock.MPS.LA.square()
            .contract(superblock.MPS.LA.square().log().eval(), idx<1>({0},{0}));
    return static_cast<double>(result1(0)) ;
}

double class_observables::get_variance(){

    if (Sim2String[sim].find("TEBD") != std::string::npos ){
        Textra::Tensor<4,Scalar> theta = superblock.MPS.get_theta();
        Textra::Tensor<4,Scalar> H_sq  = superblock.H.asTensor.contract(superblock.H.asTensor, idx<2>({0,1},{2,3}));
        Textra::Tensor<0,Scalar> Var   =  H_sq.contract(theta.conjugate(), idx<2>({0,1},{0,1}))
                .contract(theta,             idx<4>({0,1,2,3},{0,1,2,3}));
        return (static_cast<double>(Var(0)) - std::pow(get_energy(), 2));
    }else{
        return get_expectationvalue(superblock.H.M.square()) - std::pow(get_expectationvalue(superblock.H.M), 2);
    }
}

double class_observables::get_truncation_error(){
    return superblock.truncation_error;
}

long class_observables::get_chi(){
    return superblock.chi;
}

long class_observables::get_chi_max(){
    return superblock.chi_max;
}

long class_observables::get_chain_length(){
    return superblock.chain_length;
}

void class_observables::print_status_full(int verbosity) {
    if (verbosity >= 1) {
        cout << setprecision(16) << fixed << left;
        cout << setw(20) << "E_inf_" + Sim2String.at(sim) << " = "   << get_energy() << '\n';
        if(verbosity >= 2) {
            cout  << setw(20) << "Entanglement Entropy" << " = " << setprecision(16) << fixed      << get_entropy() << '\n';
            cout  << setw(20) << "Truncation error"     << " = " << setprecision(4)  << scientific << get_truncation_error() << '\n';
            cout  << setw(20) << "Variance"             << " = " << setprecision(16) << scientific << get_variance() << '\n';
            cout  << setw(20) << "chi_max"              << " = " << setprecision(4)  << fixed      << get_chi_max() << '\n';
            cout  << setw(20) << "chi    "              << " = " << setprecision(4)  << fixed      << get_chi()     << '\n';

        }
    }
}


void class_observables::print_status_update(int step) {
    cout << left  << Sim2String.at(sim)
         << left  << setw(10) << step
         << left  << "E: "            << setw(25) << setprecision(16)    << fixed   << get_energy()
         << left  << "Var(E): "       << setw(25) << setprecision(16)    << fixed   << get_variance()
         << left  << "Entanglement Entropy: "   << setw(25) << setprecision(16)    << fixed   << get_entropy()
         << left  << "SVD truncation: "   << setw(25) << setprecision(16)    << fixed   << get_truncation_error()
         << left  << "\u03C7_max: "  << setw(3)  << setprecision(4)     << fixed   << get_chi_max()
         << left  << "\u03C7: "      << setw(20) << setprecision(4)     << fixed   << get_chi()
         << '\n';
}
