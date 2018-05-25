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
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_finite_chain_sweeper.h>

using namespace std;
using namespace Textra;




class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_)
        :superblock(std::move(superblock_)), sim(sim_)
{

    mom_vec0 = superblock->H->compute_G(a0,4);
    mom_vecA = superblock->H->compute_G(a,4);
    mom_vecB = superblock->H->compute_G(b,4);
    mom_vecC = superblock->H->compute_G(c,4);
    mom_vecD = superblock->H->compute_G(d,4);
    set_profiling_labels();

}

class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_,
                                     std::shared_ptr<class_finite_chain_sweeper> env_storage_,
                                     SimulationType sim_)
        :class_measurement(superblock_,sim_)

{
        env_storage = std::move(env_storage_);
}


void class_measurement::compute_all_observables_from_superblock(){
    if(!is_measured) {
//        superblock->set_current_dimensions();
        if (superblock->chiL != superblock->chiR) { return;}
        if (superblock->MPS->LA.size() != superblock->MPS->LB.size()) { return ;}
        if (superblock->chi <= 2) { return; }

        superblock->MPS->compute_mps_components();
        energy1                       = compute_energy_MPO();
        energy2                       = compute_energy_H();
        variance1                     = compute_energy_variance_MPO();
        variance2                     = compute_energy_variance_H();
        std::tie(energy3, variance3)  = compute_infinite_moments_G(a, mom_vecA);
        std::tie(energy4, variance4)  = compute_infinite_moments_G(b, mom_vecB);
        std::tie(energy5, variance5)  = compute_infinite_moments_G(c, mom_vecC);
        std::tie(energy6, variance6)  = compute_infinite_moments_G(d, mom_vecD);
        entanglement_entropy          = compute_entanglement_entropy();
        truncation_error              = superblock->MPS->truncation_error;
//        parity                        = compute_parity();
        is_measured = true;
    }
}





double class_measurement::compute_energy_MPO(){
    t_ene_mpo.tic();

    if(sim == SimulationType::iDMRG ){
        Tensor<Scalar, 0>  E_two_sites =
                superblock->Lblock->block
                        .contract(superblock->MPS->theta,                     idx({0},{1}))
                        .contract(superblock->HA->MPO,                        idx({1,2},{0,2}))
                        .contract(superblock->HB->MPO,                        idx({3,1},{0,2}))
                        .contract(superblock->MPS->theta.conjugate(),         idx({0,2,4},{1,0,2}))
                        .contract(superblock->Rblock->block,                  idx({0,2,1},{0,1,2}));
        t_ene_mpo.toc();
        return std::real(E_two_sites(0) / 2.0 );

    }else{
        Tensor<Scalar, 0>  E_all_sites =
                superblock->Lblock->block
                        .contract(superblock->MPS->theta,                     idx({0},{1}))
                        .contract(superblock->HA->MPO,                        idx({1,2},{0,2}))
                        .contract(superblock->HB->MPO,                        idx({3,1},{0,2}))
                        .contract(superblock->MPS->theta.conjugate(),         idx({0,2,4},{1,0,2}))
                        .contract(superblock->Rblock->block,                  idx({0,2,1},{0,1,2}));
        double L = superblock->Lblock->size + superblock->Rblock->size + 2;
        t_ene_mpo.toc();
        energy_all_sites = std::real(E_all_sites(0));
        return std::real(energy_all_sites)/L;
    }
}


double class_measurement::compute_energy_H(){
    t_ene_ham.tic();
    E_evn = superblock->MPS->theta_evn_normalized
                .contract(Matrix_to_Tensor(superblock->H->h[0],2,2,2,2),     idx({0, 2}, {0, 1}))
                .contract(superblock->MPS->theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
                .contract(superblock->MPS->l_evn,                            idx({0, 2}, {0, 1}))
                .contract(superblock->MPS->r_evn,                            idx({0, 1}, {0, 1}));
    E_odd  = superblock->MPS-> theta_odd_normalized
                    .contract(Matrix_to_Tensor(superblock->H->h[1],2,2,2,2),    idx({0, 2}, {0, 1}))
                    .contract(superblock->MPS->theta_odd_normalized.conjugate(),idx({2, 3}, {0,2}))
                    .contract(superblock->MPS->l_odd,                           idx({0,2}, {0,1}))
                    .contract(superblock->MPS->r_odd,                           idx({0,1},{0,1}));
    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );
    t_ene_ham.toc();
    return 0.5*std::real(E_evn(0) + E_odd(0)) ;
}



double class_measurement::compute_entanglement_entropy(){
    t_entropy.tic();
    Tensor<Scalar,0> SA  = -superblock->MPS->LA.square()
            .contract(superblock->MPS->LA.square().log().eval(), idx({0},{0}));
    t_entropy.toc();
    return std::real(SA(0));
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


double class_measurement::compute_finite_energy(){
    Eigen::Tensor<Scalar,3> L = env_storage->get_ENV_L().front().block;
    Eigen::TensorRef<Eigen::Tensor<Scalar,3>> temp;

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto mpoL  = env_storage->get_MPO_L().begin();
    auto endL  = env_storage->get_MPS_L().end();
    int iter = 0;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1> &LB_left = std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3> &GA      = std::get<1>(*mpsL);
        assert(LB_left.dimension(0) == L.dimension(0));
        assert(LB_left.dimension(0) == GA.dimension(1));

        temp = L.contract(asDiagonal(LB_left), idx({0},{0}))
                .contract(asDiagonal(LB_left), idx({0},{0}))
                .contract(mpoL->MPO,           idx({0},{0}))
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
        const Eigen::Tensor<Scalar,3> &GB      = std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB      = std::get<1>(*mpsR);
        assert(GB.dimension(1) == L.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));
        temp = L.contract(GB,            idx({0},{1}))
                .contract(GB.conjugate(),idx({0},{1}))
                .contract(mpoR->MPO,     idx({0,1,3},{0,2,3}))
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
    return std::real(E_all_sites(0));
}


double class_measurement::compute_finite_norm(){

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto endL  = env_storage->get_MPS_L().end();
    const Eigen::Tensor<Scalar,1>  & LB_left = std::get<0>(*mpsL);
    const Eigen::Tensor<Scalar,3>  & GA      = std::get<1>(*mpsL);
    Eigen::Tensor<Scalar,2> chain =
             asDiagonal(LB_left)
            .contract(asDiagonal(LB_left), idx({0},{1}))
            .contract(GA                 , idx({1},{1}))
            .contract(GA.conjugate()     , idx({0,1},{1,0}));
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;
    mpsL++;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1>  & LB_left = std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3>  & GA      = std::get<1>(*mpsL);
        assert(LB_left.dimension(0) == chain.dimension(0));
        assert(LB_left.dimension(0) == GA.dimension(1));
        temp = chain
                .contract(asDiagonal(LB_left), idx({0},{0}))
                .contract(asDiagonal(LB_left), idx({0},{0}))
                .contract(GA                 , idx({0},{1}))
                .contract(GA.conjugate()     , idx({0,1},{1,0}));
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
        const Eigen::Tensor<Scalar,3> &GB  = std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB  = std::get<1>(*mpsR);
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
    Scalar norm = Textra::Tensor2_to_Matrix(chain).trace();
    std::cout << setprecision(16) << "Norm: " << norm << std::endl;
    assert(std::abs(norm - 1.0) < 1e-10 );
    return std::abs(norm);
}


void class_measurement::compute_finite_mps_state(){

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto endL  = env_storage->get_MPS_L().end();
    Eigen::Tensor<Scalar,2> chain(1,1);
    chain.setConstant(1.0);
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;

    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1>  & LB_left = std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3>  & GA      = std::get<1>(*mpsL);
        assert(LB_left.dimension(0) == GA.dimension(1));

        temp = chain
                .contract(asDiagonal(LB_left), idx({1},{0}))
                .contract(GA                 , idx({1},{1}))
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
        const Eigen::Tensor<Scalar,3> &GB  = std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB  = std::get<1>(*mpsR);
        assert(LB.dimension(0) == GB.dimension(2));
        temp = chain
                .contract(GB              , idx({1},{1}))
                .contract(asDiagonal(LB)  , idx({2},{0}))
                .reshape(array2{GB.dimension(0) * chain.dimension(0), GB.dimension(2)});
        chain = temp;
        mpsR++;
    }
    Scalar norm = Textra::Tensor2_to_Matrix(chain).norm();
    std::cout << "The state: \n" << chain << std::endl;
    std::cout << setprecision(16) << "Norm of state: " << norm << std::endl;
    assert(std::abs(norm - 1.0) < 1e-10 );
//    return std::abs(norm);
}


double class_measurement::get_energy1(){return energy1;};
double class_measurement::get_energy2(){return energy2;};
double class_measurement::get_energy3(){return energy3;};
double class_measurement::get_energy4(){return energy4;};
double class_measurement::get_energy5(){return energy5;};
double class_measurement::get_energy6(){return energy6;};

double class_measurement::get_variance1(){return variance1;};
double class_measurement::get_variance2(){return variance2;};
double class_measurement::get_variance3(){return variance3;};
double class_measurement::get_variance4(){return variance4;};
double class_measurement::get_variance5(){return variance5;};
double class_measurement::get_variance6(){return variance6;};

double class_measurement::get_entanglement_entropy(){
    return entanglement_entropy;
}

double class_measurement::get_truncation_error(){
    return truncation_error;
}

double class_measurement::get_parity(){
    return parity;
}


long class_measurement::get_chi(){
    return superblock->chi;
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

