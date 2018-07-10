//
// Created by david on 2017-11-12.
//



#include <iomanip>
#include <complex>
#include <mps_routines/class_measurement.h>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_hamiltonian.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <general/nmspc_quantum_mechanics.h>

using namespace std;
using namespace Textra;




class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_)
        :superblock(std::move(superblock_)), sim_type(sim_)
{

    set_profiling_labels();
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    h_evn = superblock->HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    h_odd = superblock->HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    mom_vecA = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);

}

class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_,
                                     std::shared_ptr<class_finite_chain_sweeper> env_storage_,
                                     SimulationType sim_)
        :class_measurement(superblock_,sim_)

{
        env_storage = std::move(env_storage_);
}


void class_measurement::compute_all_observables_from_superblock(){
    if (is_measured){return;}
    if (superblock->MPS->chiA() != superblock->MPS->chiB()) { return;}
    if (superblock->MPS->chiA() != superblock->MPS->chiC()) { return ;}
    if (superblock->MPS->chiC() <= 2) { return; }

    switch(sim_type){
        case SimulationType::iDMRG:
            superblock->MPS->compute_mps_components();
            compute_energy_mpo();
            compute_energy_variance_mpo();
            compute_energy_ham();
            compute_energy_variance_ham();
            compute_energy_and_variance_mom(a, mom_vecA);
            break;
        case SimulationType::xDMRG:
        case SimulationType::fDMRG:
            superblock->MPS->theta = superblock->MPS->get_theta();
            compute_energy_mpo();
            compute_energy_variance_mpo();
            break;
        case SimulationType::iTEBD:
            superblock->MPS->compute_mps_components();
            compute_energy_ham();
            compute_energy_variance_ham();
            compute_energy_and_variance_mom(a, mom_vecA);
            break;
    }

    compute_entanglement_entropy();
//    compute_parity();
    is_measured = true;

}


void class_measurement::compute_all_observables_from_finite_chain(){
    compute_finite_chain_energy();
    compute_finite_chain_energy_variance();
//    compute_finite_chain_mps_state(); // Runs out of memory for L > 20 or so, because the wave vector is doubled in size for each additional site.
    compute_finite_chain_norm();
}




void class_measurement::compute_energy_mpo(){
    t_ene_mpo.tic();
    double         L = superblock->Lblock->size + superblock->Rblock->size + 2;
    Eigen::Tensor<Scalar, 0>  E =
            superblock->Lblock->block
                    .contract(superblock->MPS->theta,                     idx({0},{1}))
                    .contract(superblock->HA->MPO,                        idx({1,2},{0,2}))
                    .contract(superblock->HB->MPO,                        idx({3,1},{0,2}))
                    .contract(superblock->MPS->theta.conjugate(),         idx({0,2,4},{1,0,2}))
                    .contract(superblock->Rblock->block,                  idx({0,2,1},{0,1,2}));
    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");
    if(sim_type == SimulationType::iDMRG ){
        //iDMRG uses reduced mpo in the environments!
        energy_mpo_all_sites = std::real(E(0))/2.0*L;
        energy_mpo           = std::real(E(0))/2.0;

    }else{
        energy_mpo_all_sites = std::real(E(0));
        energy_mpo           = std::real(energy_mpo_all_sites)/L;
    }

    t_ene_mpo.toc();

}


void class_measurement::compute_energy_ham(){
    t_ene_ham.tic();
    E_evn = superblock->MPS->theta_evn_normalized
                .contract(Matrix_to_Tensor(h_evn,2,2,2,2),     idx({0, 2}, {0, 1}))
                .contract(superblock->MPS->theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
                .contract(superblock->MPS->l_evn,                            idx({0, 2}, {0, 1}))
                .contract(superblock->MPS->r_evn,                            idx({0, 1}, {0, 1}));
    E_odd  = superblock->MPS-> theta_odd_normalized
                    .contract(Matrix_to_Tensor(h_odd,2,2,2,2),    idx({0, 2}, {0, 1}))
                    .contract(superblock->MPS->theta_odd_normalized.conjugate(),idx({2, 3}, {0,2}))
                    .contract(superblock->MPS->l_odd,                           idx({0,2}, {0,1}))
                    .contract(superblock->MPS->r_odd,                           idx({0,1},{0,1}));
    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );

    t_ene_ham.toc();
    energy_ham = 0.5*std::real(E_evn(0) + E_odd(0));
}



void class_measurement::compute_entanglement_entropy(){
    t_entropy.tic();
    Eigen::Tensor<Scalar,0> SA  = -superblock->MPS->LC.square()
            .contract(superblock->MPS->LC.square().log().eval(), idx({0},{0}));
    entanglement_entropy =  std::real(SA(0));
    t_entropy.toc();
}

//
//double class_measurement::compute_parity(){
////    t_ene_mpo.tic();
//    Tensor<T, 0>  E_two_site = superblock->MPS->theta
//                    .contract(parity_mpo.MPO,     idx({0},{0}))
//                    .contract(parity_mpo.MPO,     idx({1},{0}))
//                    .contract(superblock->MPS->theta.conjugate(),         idx({0,1,2,3},{1,3,0,2}));
////    t_ene_mpo.toc();
//    return std::real(E_two_site(0) / 2.0 );
//}


void class_measurement::compute_finite_chain_energy(){
    Eigen::Tensor<Scalar,3> L = env_storage->get_ENV_L().front().block;
    Eigen::TensorRef<Eigen::Tensor<Scalar,3>> temp;

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto mpoL  = env_storage->get_MPO_L().begin();
    auto endL  = env_storage->get_MPS_L().end();
    int iter = 0;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1> &LA = mpsL->get_L(); // std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3> &GA = mpsL->get_G(); // std::get<1>(*mpsL);
        assert(LA.dimension(0) == L.dimension(0));
        assert(LA.dimension(0) == GA.dimension(1));

        temp = L.contract(asDiagonal(LA), idx({0},{0}))
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(mpoL->get()->MPO,idx({0},{0}))
                .contract(GA,                  idx({0,3},{1,0}))
                .contract(GA.conjugate(),      idx({0,2},{1,0}))
                .shuffle(array3{1,2,0});
        L = temp;
        mpsL++;
        mpoL++;
        iter++;
    }


    //Contract the center point
    auto &MPS_C = env_storage->get_MPS_C();
    temp = L.contract(asDiagonal(MPS_C) , idx({0},{0}))
            .contract(asDiagonal(MPS_C) , idx({0},{0}))
            .shuffle(array3{1,2,0});
    L = temp;

    //Contract the right half of the chain
    Eigen::Tensor<Scalar,3> R = env_storage->get_ENV_R().back().block;
    auto mpsR  = env_storage->get_MPS_R().begin();
    auto mpoR  = env_storage->get_MPO_R().begin();
    auto endR  = env_storage->get_MPS_R().end();
    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB = mpsR->get_G(); // std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB = mpsR->get_L(); // std::get<1>(*mpsR);
        assert(GB.dimension(1) == L.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));
        temp = L.contract(GB,            idx({0},{1}))
                .contract(GB.conjugate(),idx({0},{1}))
                .contract(mpoR->get()->MPO,     idx({0,1,3},{0,2,3}))
                .contract(asDiagonal(LB),idx({0},{0}))
                .contract(asDiagonal(LB),idx({0},{0}))
                .shuffle(array3{1,2,0});
        L = temp;
        mpsR++;
        mpoR++;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar,0> E_all_sites = L.contract(R, idx({0,1,2},{0,1,2}));
    energy_chain = std::real(E_all_sites(0));
    std::cout << setprecision(16) << "E all sites: " << energy_chain << " | per site: " << energy_chain / env_storage->get_length() << std::endl;
}


void class_measurement::compute_finite_chain_norm(){

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto endL  = env_storage->get_MPS_L().end();
    const Eigen::Tensor<Scalar,1>  & LA_ = mpsL->get_L(); // std::get<0>(*mpsL);
    const Eigen::Tensor<Scalar,3>  & GA_ = mpsL->get_G(); // std::get<1>(*mpsL);
    Eigen::Tensor<Scalar,2> chain =
             asDiagonal(LA_)
            .contract(asDiagonal(LA_), idx({0},{1}))
            .contract(GA_                 , idx({1},{1}))
            .contract(GA_.conjugate()     , idx({0,1},{1,0}));
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;
    mpsL++;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1>  & LA = mpsL->get_L() ; // std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3>  & GA = mpsL->get_G() ; // std::get<1>(*mpsL);
        assert(LA.dimension(0) == chain.dimension(0));
        assert(LA.dimension(0) == GA.dimension(1));
        temp = chain
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(GA            , idx({0},{1}))
                .contract(GA.conjugate(), idx({0,1},{1,0}));
        chain = temp;
        mpsL++;
    }

//    Contract the center point
    auto &MPS_C = env_storage->get_MPS_C();
    temp = chain
            .contract(asDiagonal(MPS_C), idx({0},{0}))
            .contract(asDiagonal(MPS_C), idx({0},{0}));
    chain = temp;

    //Contract the right half of the chain
    auto mpsR  = env_storage->get_MPS_R().begin();
    auto endR  = env_storage->get_MPS_R().end();

    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB  = mpsR->get_G(); //std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB  = mpsR->get_L(); //std::get<1>(*mpsR);
        assert(GB.dimension(1) == chain.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));

        temp = chain
                .contract(GB              , idx({0},{1}))
                .contract(GB.conjugate()  , idx({0,1},{1,0}))
                .contract(asDiagonal(LB)  , idx({0},{0}))
                .contract(asDiagonal(LB)  , idx({0},{0}));
        chain = temp;
        mpsR++;
    }
    norm_chain = std::abs(Textra::Tensor2_to_Matrix(chain).trace());
    std::cout << setprecision(16) << "Norm: " << norm_chain << std::endl;
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
}


void class_measurement::compute_finite_chain_mps_state(){

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto endL  = env_storage->get_MPS_L().end();
    Eigen::Tensor<Scalar,2> chain(1,1);
    chain.setConstant(1.0);
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;

    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1>  & LA = mpsL->get_L(); //std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3>  & GA = mpsL->get_G(); //std::get<1>(*mpsL);
        assert(LA.dimension(0) == GA.dimension(1));

        temp = chain
                .contract(asDiagonal(LA), idx({1},{0}))
                .contract(GA            , idx({1},{1}))
                .reshape(array2{GA.dimension(0) * chain.dimension(0), GA.dimension(2)});
        chain = temp;
        mpsL++;
    }

//    Contract the center point
    auto &MPS_C = env_storage->get_MPS_C();
    temp = chain.contract(asDiagonal(MPS_C), idx({1},{0}));
    chain = temp;

    //Contract the right half of the chain
    auto mpsR  = env_storage->get_MPS_R().begin();
    auto endR  = env_storage->get_MPS_R().end();

    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB  = mpsR->get_G(); // std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB  = mpsR->get_L(); // std::get<1>(*mpsR);
        assert(LB.dimension(0) == GB.dimension(2));
        temp = chain
                .contract(GB              , idx({1},{1}))
                .contract(asDiagonal(LB)  , idx({2},{0}))
                .reshape(array2{GB.dimension(0) * chain.dimension(0), GB.dimension(2)});
        chain = temp;
        mpsR++;
    }
    mps_chain = chain.reshape(array1{chain.dimension(0)});
    norm_chain = Textra::Tensor2_to_Matrix(chain).norm();
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
}


double class_measurement::get_energy_mpo(){return energy_mpo;};
double class_measurement::get_energy_ham(){return energy_ham;};
double class_measurement::get_energy_mom(){return energy_mom;};

double class_measurement::get_variance_mpo(){return variance_mpo;};
double class_measurement::get_variance_ham(){return variance_ham;};
double class_measurement::get_variance_mom(){return variance_mom;};

double class_measurement::get_entanglement_entropy(){
    compute_entanglement_entropy();
    return entanglement_entropy;
}

double class_measurement::get_truncation_error(){
    return superblock->MPS->truncation_error;
}

double class_measurement::get_parity(){
    return parity;
}


long class_measurement::get_chi(){
    return superblock->MPS->chiC();
}

long class_measurement::get_chain_length(){
    return superblock->environment_size + 2;
}


// Profiling

void class_measurement::set_profiling_labels() {
    using namespace settings::profiling;
    t_ene_mpo.set_properties(on, precision,"↳ Energy (MPO)           ");
    t_ene_ham.set_properties(on, precision,"↳ Energy (Ham)           ");
    t_ene_gen.set_properties(on, precision,"↳ Energy (Gen)           ");
    t_var_mpo.set_properties(on, precision,"↳ Variance (MPO)         ");
    t_var_ham.set_properties(on, precision,"↳ Variance (Ham)         ");
    t_var_gen.set_properties(on, precision,"↳ Variance (Gen)         ");
    t_entropy.set_properties(on, precision,"↳ Ent. Entropy           ");
    t_temp1.set_properties(on, precision,  "↳ Temp1                  ");
    t_temp2.set_properties(on, precision,  "↳ Temp2                  ");
    t_temp3.set_properties(on, precision,  "↳ Temp3                  ");
    t_temp4.set_properties(on, precision,  "↳ Temp4                  ");

}

void class_measurement::print_profiling(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\nComputing observables breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_ene_mpo.print_time_w_percent(t_parent);
        t_ene_ham.print_time_w_percent(t_parent);
        t_ene_gen.print_time_w_percent(t_parent);
        t_var_mpo.print_time_w_percent(t_parent);
        t_var_ham.print_time_w_percent(t_parent);
        t_var_gen.print_time_w_percent(t_parent);
        t_entropy.print_time_w_percent(t_parent);
        t_temp1.print_time_w_percent(t_parent);
        t_temp2.print_time_w_percent(t_parent);
        t_temp3.print_time_w_percent(t_parent);
        t_temp4.print_time_w_percent(t_parent);
    }
}

