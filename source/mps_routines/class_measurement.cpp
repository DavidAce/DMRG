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
#include <general/class_arpackpp_wrapper.h>
#include <algorithms/class_base_algorithm.h>

using namespace std;
using namespace Textra;




class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_)
        :superblock(std::move(superblock_)), sim(sim_)
{

}


double class_measurement::first_moment(){
    superblock->update_bond_dimensions();
    if (superblock->chiL != superblock->chiR) { return 1.0;}
    if (superblock->MPS->LA.size() != superblock->MPS->LB.size()) { return 1.0;}
    if (superblock->chi <= 2) { return 1.0;}
    using T = std::complex<double>;
    class_arpackpp_wrapper eig;
    Textra::array<4> shape4 = superblock->shape4;
    long sizeL     = shape4[1] * shape4[1];
    long sizeR     = shape4[3] * shape4[3];
    Textra::Tensor<T,4> thetaR      = superblock->MPS->thetaR().cast<T>();
    Textra::Tensor<T,3> A           = superblock->MPS->A().cast<T>();
    Textra::Tensor<T,2> transf_G    = thetaR.contract(superblock->H->G, Textra::idx<2>({0,1},{0,1}))
                                       .contract(thetaR.conjugate(), Textra::idx<2>({2,3},{0,1}))
                                       .shuffle(Textra::array4{0,2,1,3})
                                       .reshape(Textra::array2{sizeL,sizeR}) ;
//    std::cout << "dimensions:" << A.dimensions() << std::endl;
//    assert(A.dimension(1) == A.dimension(2) and "Dimensions mismatch of A matrix");
    Textra::Tensor<T,2> transf_ID = A.contract(A.conjugate(), Textra::idx<1>({0},{0}))
        .shuffle(Textra::array4{0,2,1,3})
        .reshape(Textra::array2{sizeL,sizeR}) ;

    if(eigvecG.size() == transf_ID.dimension(0)) {
        std::tie(eigvecG, lambdaG) = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, true>(transf_G.data() , sizeL,sizeR, 1, eigvecG.data());
        std::tie(eigvecID,lambdaID)= eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, true>(transf_ID.data(), sizeL,sizeR, 1, eigvecID.data());
    }else{
        std::tie(eigvecG, lambdaG) = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, true>(transf_G.data() , sizeL,sizeR, 1);
        std::tie(eigvecID,lambdaID)= eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, true>(transf_ID.data(), sizeL,sizeR, 1);
    }

    T G             = (lambdaG(0));
    T ID            = (lambdaID(0));
    T G_ID          = G/ID;   //
    T sqrtG_ID      = std::sqrt(G)/ID; //Gives spikes
    T sqrtG_sqrtID  = std::sqrt(G_ID);// --> best yet. It is very stable and doesnt give spikes. It should be squared later, as in the paper


    switch(sim){
        case SimulationType::iTEBD:
        case SimulationType ::FES_iTEBD:
            variance1 = 1;
            break;
        default:
            variance1 = get_full_variance();
    }

    variance2 = std::abs(std::log(std::pow(std::abs(sqrtG_sqrtID),2)));
    variance3 = variance2;//std::abs(std::log(Fa2) + std::log(Fa1));

    //NOTES:
    //  std::abs(std::log(Fa2) + std::log(Fa1)); // gives exactly the same result as      std::abs(std::log(std::pow(std::abs(sqrtG_sqrtID),2)));


    return variance1;
}


double class_measurement::get_full_variance(){
    Tensor<double,4> theta  = superblock->MPS->get_theta();
    Tensor<double,0> E =
            superblock->Lblock->block
                    .contract(theta.conjugate(),        idx<1>({0},{2}))
                    .contract(superblock->H->MM ,       idx<3>({1,2,3},{0,2,3}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(theta,                    idx<3>({0,3,4},{2,0,1}))
                    .contract(superblock->Rblock->block,idx<3>({0,2,1},{0,1,2}));

    Tensor<double,0> E2 =
            superblock->Lblock2->block
                    .contract(theta.conjugate()  ,         idx<1>({0},{2}))
                    .contract(superblock->H->MMMM,         idx<4>({2,1,3,4},{2,0,4,5}))
                    .contract(theta,                       idx<3>({0,4,5},{2,0,1}))
                    .contract(superblock->Rblock2->block,  idx<4>({0,3,1,2},{0,1,2,3}));

    return std::abs(E2(0) - E(0)*E(0))/ superblock->chain_length / superblock->chain_length;

}

double class_measurement::get_expectationvalue(const Tensor<double,4> &MPO){
    Tensor<double,4> theta  = superblock->MPS->get_theta();
    Tensor<double,0> result =
            superblock->Lblock->block
                    .contract(theta.conjugate(),idx<1>({0},{2}))
                    .contract(MPO,              idx<2>({1,2},{0,2}))
                    .contract(MPO,              idx<2>({3,1},{0,2}))
                    .contract(theta,            idx<3>({0,2,4},{2,0,1}))
                    .contract(superblock->Rblock->block,     idx<3>({0,1,2},{0,2,1}));
    return result(0)/ superblock->chain_length;
}

double class_measurement::get_expectationvalue_sq(const Tensor<double,4> &MPO){
    Tensor<double,4> theta  = superblock->MPS->get_theta();
    Tensor<double,0> result =
            superblock->Lblock2->block
                    .contract(theta.conjugate(),idx<1>({0}  ,{2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,2},{0,2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,3},{0,2}))
                    .contract(theta,            idx<3>({0,3,5},{2,0,1}))
                    .contract(superblock->Rblock2->block,     idx<4>({0,1,2,3},{0,2,3,1}));
    return result(0) / superblock->chain_length / superblock->chain_length;
}

double class_measurement::get_expectationvalue(const Tensor<std::complex<double>,4> &MPO){
    using T = std::complex<double>;
    Tensor<T,4> theta  = superblock->MPS->get_theta().cast<T>();
    Tensor<T,0> result =
            superblock->Lblock->block.cast<T>()
                    .contract(theta.conjugate(),idx<1>({0},{2}))
                    .contract(MPO,              idx<2>({1,2},{0,2}))
                    .contract(MPO,              idx<2>({3,1},{0,2}))
                    .contract(theta,            idx<3>({0,2,4},{2,0,1}))
                    .contract(superblock->Rblock->block.cast<T>(),     idx<3>({0,1,2},{0,2,1}));
    return result(0).real()/ superblock->chain_length;
}

double class_measurement::get_expectationvalue_sq(const Tensor<std::complex<double>,4> &MPO){
    using T = std::complex<double>;
    Tensor<T,4> theta  = superblock->MPS->get_theta().cast<T>();
    Tensor<T,0> result =
            superblock->Lblock2->block.cast<T>()
                    .contract(theta.conjugate(),idx<1>({0}  ,{2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,2},{0,2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,3},{0,2}))
                    .contract(theta,            idx<3>({0,3,5},{2,0,1}))
                    .contract(superblock->Rblock2->block.cast<T>(),     idx<4>({0,1,2,3},{0,2,3,1}));
    return result(0).real() / superblock->chain_length / superblock->chain_length;
}


double class_measurement::get_second_cumulant(){
    return first_moment();
}



double class_measurement::get_energy(){
    Tensor<double,4> theta  = superblock->MPS->get_theta();
    Tensor<double,0> EA,EB;
    switch (sim){
        case SimulationType::iDMRG:
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
        case SimulationType::FES_iDMRG:
            EA =
                    superblock->Lblock->block
                            .contract(theta.conjugate(),        idx<1>({0},{2}))
                            .contract(superblock->H->MM ,       idx<3>({1,2,3},{0,2,3}))//  idx<3>({1,2,3},{0,4,5}))
                            .contract(theta,                    idx<3>({0,3,4},{2,0,1}))
                            .contract(superblock->Rblock->block,idx<3>({0,2,1},{0,1,2}));
//            theta = superblock->MPS->get_theta_swapped();
//            EB =
//                    superblock->Lblock->block
//                            .contract(theta.conjugate(),        idx<1>({0},{2}))
//                            .contract(superblock->H->MM ,       idx<3>({1,2,3},{0,2,3}))//  idx<3>({1,2,3},{0,4,5}))
//                            .contract(theta,                    idx<3>({0,3,4},{2,0,1}))
//                            .contract(superblock->Rblock->block,idx<3>({0,2,1},{0,1,2}));
//        std::cout << setprecision(8);
//            std::cout << "\nEA:" << EA(0)/ superblock->chain_length << std::endl;
//            std::cout << "EB:" << EB(0)/ superblock->chain_length << std::endl;


//            return 0.5*(EA(0)+EB(0))/ superblock->chain_length;
            return EA(0)/ superblock->chain_length;

        case SimulationType::iTEBD:
        case SimulationType::FES_iTEBD:
            EA = superblock->H->H_asTensor
                    .contract(theta.conjugate(), idx<2>({0, 1}, {0, 1}))
                    .contract(theta, idx<4>({0, 1, 2, 3}, {0, 1, 2, 3}));

            theta = superblock->MPS->get_theta_swapped();
            EB = superblock->H->H_asTensor
                    .contract(theta.conjugate(), idx<2>({0, 1}, {0, 1}))
                    .contract(theta, idx<4>({0, 1, 2, 3}, {0, 1, 2, 3}));
            return static_cast<double>(0.5 * (EA(0) + EB(0))) ;

    }
}

double class_measurement::get_entanglement_entropy(){
    Tensor<Scalar,0> SA  = -superblock->MPS->LA.square()
            .contract(superblock->MPS->LA.square().log().eval(), idx<1>({0},{0}));
    Tensor<Scalar,0> SB  = -superblock->MPS->LB.square()
            .contract(superblock->MPS->LB.square().log().eval(), idx<1>({0},{0}));
    return static_cast<double>(0.5*(SA(0)+SB(0))) ;
}


double class_measurement::get_variance(){
    return get_second_cumulant();
}

double class_measurement::get_variance1(){return variance1;};
double class_measurement::get_variance2(){return variance2;};
double class_measurement::get_variance3(){return variance3;};

double class_measurement::get_truncation_error(){
    return superblock->truncation_error;
}

long class_measurement::get_chi(){
    return superblock->chi;
}

long class_measurement::get_chain_length(){
    return superblock->chain_length;
}

