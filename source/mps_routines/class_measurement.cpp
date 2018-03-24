//
// Created by david on 2017-11-12.
//



#include <iomanip>
#include <complex>
#include <mps_routines/class_measurement.h>
#include <general/nmspc_tensor_extra.h>
#include <Eigen/Eigenvalues>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <algorithms/class_base_algorithm.h>

using namespace std;
using namespace Textra;




class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_)
        :superblock(std::move(superblock_)), sim(sim_)
{

    mom_vecA = superblock->H->compute_G(a,4);
    mom_vecB = superblock->H->compute_G(b,4);
    mom_vecC = superblock->H->compute_G(c,4);
    mom_vecD = superblock->H->compute_G(d,4);

}


void class_measurement::do_full_measurement(){
    if(!is_measured) {
        superblock->set_current_dimensions();
        if (superblock->chiL != superblock->chiR) { return;}
        if (superblock->MPS->LA.size() != superblock->MPS->LB.size()) { return ;}
        if (superblock->chi <= 2) { return; }

        superblock->MPS->compute_mps_components();
        energy1                       = compute_energy_MPO();
        energy2                       = compute_energy_H();
        variance1                     = compute_infinite_variance_MPO();
        variance2                     = compute_infinite_variance_H();
        std::tie(energy3, variance3)  = compute_infinite_moments_G(a, mom_vecA);
        std::tie(energy4, variance4)  = compute_infinite_moments_G(b, mom_vecB);
        std::tie(energy5, variance5)  = compute_infinite_moments_G(c, mom_vecC);
        std::tie(energy6, variance6)  = compute_infinite_moments_G(d, mom_vecD);
        entanglement_entropy          = compute_entanglement_entropy();
        is_measured = true;
    }
}


double class_measurement::compute_energy_MPO(){
    Tensor<Scalar, 0>  E_two_site =
            superblock->Lblock->block
            .contract(superblock->MPS->theta,                     idx({0},{1}))
            .contract(superblock->H->M,                           idx({1,2},{0,2}))
            .contract(superblock->H->M,                           idx({3,1},{0,2}))
            .contract(superblock->MPS->theta.conjugate(),         idx({0,2,4},{1,0,2}))
            .contract(superblock->Rblock->block,                  idx({0,2,1},{0,1,2}));
    return std::real(E_two_site(0) / 2.0 );
}


double class_measurement::compute_energy_H(){
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
    return 0.5*std::real(E_evn(0) + E_odd(0)) ;

}



double class_measurement::compute_entanglement_entropy(){
    Tensor<Scalar,0> SA  = -superblock->MPS->LA.square()
            .contract(superblock->MPS->LA.square().log().eval(), idx({0},{0}));
    Tensor<Scalar,0> SB  = -superblock->MPS->LB.square()
            .contract(superblock->MPS->LB.square().log().eval(), idx({0},{0}));
    return std::real(0.5*(SA(0)+SB(0))) ;
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
    return superblock->MPS->truncation_error;
}

long class_measurement::get_chi(){
    return superblock->chi;
}

long class_measurement::get_chain_length(){
    return superblock->chain_length;
}

