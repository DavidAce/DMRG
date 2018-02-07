//
// Created by david on 2017-11-12.
//

#include "class_measurement.h"
#include <iomanip>
#include <complex>
#include <Eigen/Eigenvalues>
#include <mps_routines/class_superblock.h>
#include <general/class_eig_arpack_wrapper.h>


using namespace std;
using namespace Textra;




class_measurement::class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_)
        :superblock(std::move(superblock_)), sim(sim_)
{

}


double class_measurement::first_moment(){
    superblock->update_bond_dimensions();
    if (superblock->MPS.LA.size() != superblock->MPS.LB.size()) { return 1.0;}
    double h = 0.0000001;

    double a = 0.0;
    double a1  = a - h;

    using T = std::complex<double>;
    class_eig_arpack eig;

    Textra::array<4> shape4 = superblock->shape4;
    long sizeL     = shape4[1] * shape4[1];
    long sizeR     = shape4[3] * shape4[3];

    Textra::Tensor<T,4> theta = superblock->MPS.thetaR().cast<T>();
    Textra::Tensor<T,3> A     = superblock->MPS.A().cast<T>();
    Textra::Tensor<T,2> transf_Ga1 = theta.contract(superblock->H.compute_G(a1), Textra::idx<2>({0,1},{0,1})).contract(theta.conjugate(), Textra::idx<2>({2,3},{0,1})).shuffle(Textra::array4{0,2,1,3}).reshape(Textra::array2{sizeL,sizeR}) ;


    Textra::Tensor<T,2> transf_ID = A.contract(A.conjugate(), Textra::idx<1>({0},{0}))
            .shuffle(Textra::array4{0,2,1,3})
            .reshape(Textra::array2{sizeL,sizeR}) ;

    Textra::Tensor<T,1> lambda_Ga1     = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Ga1, 1);
    T ga1  = (lambda_Ga1(0));


    double H2 = get_expectationvalue_sq(superblock->H.M);
    double E2 = get_expectationvalue(superblock->H.M) * get_expectationvalue(superblock->H.M);


    variance1 = H2-E2;
    variance2 = std::fabs((std::log(std::pow(std::fabs(ga1),2))));
    return std::log(std::pow(std::abs(ga1),2));
}


double class_measurement::get_expectationvalue(const Tensor<double,4> &MPO){
    Tensor<double,4> theta  = superblock->MPS.get_theta();
    Tensor<double,0> result =
            superblock->Lblock.block
                    .contract(theta.conjugate(),idx<1>({0},{2}))
                    .contract(MPO,              idx<2>({1,2},{0,2}))
                    .contract(MPO,              idx<2>({3,1},{0,2}))
                    .contract(theta,            idx<3>({0,2,4},{2,0,1}))
                    .contract(superblock->Rblock.block,     idx<3>({0,1,2},{0,2,1}));
    return result(0)/ superblock->chain_length;
}

double class_measurement::get_expectationvalue_sq(const Tensor<double,4> &MPO){
    Tensor<double,4> theta  = superblock->MPS.get_theta();
    Tensor<double,0> result =
            superblock->Lblock2.block
                    .contract(theta.conjugate(),idx<1>({0}  ,{2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,2},{0,2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,3},{0,2}))
                    .contract(theta,            idx<3>({0,3,5},{2,0,1}))
                    .contract(superblock->Rblock2.block,     idx<4>({0,1,2,3},{0,2,3,1}));
    return result(0) / superblock->chain_length / superblock->chain_length;
}

double class_measurement::get_expectationvalue(const Tensor<std::complex<double>,4> &MPO){
    using T = std::complex<double>;
    Tensor<T,4> theta  = superblock->MPS.get_theta().cast<T>();
    Tensor<T,0> result =
            superblock->Lblock.block.cast<T>()
                    .contract(theta.conjugate(),idx<1>({0},{2}))
                    .contract(MPO,              idx<2>({1,2},{0,2}))
                    .contract(MPO,              idx<2>({3,1},{0,2}))
                    .contract(theta,            idx<3>({0,2,4},{2,0,1}))
                    .contract(superblock->Rblock.block.cast<T>(),     idx<3>({0,1,2},{0,2,1}));
    return result(0).real()/ superblock->chain_length;
}

double class_measurement::get_expectationvalue_sq(const Tensor<std::complex<double>,4> &MPO){
    using T = std::complex<double>;
    Tensor<T,4> theta  = superblock->MPS.get_theta().cast<T>();
    Tensor<T,0> result =
            superblock->Lblock2.block.cast<T>()
                    .contract(theta.conjugate(),idx<1>({0}  ,{2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,2},{0,2}))
                    .contract(MPO,              idx<2>({1,3},{0,2}))
                    .contract(MPO,              idx<2>({4,3},{0,2}))
                    .contract(theta,            idx<3>({0,3,5},{2,0,1}))
                    .contract(superblock->Rblock2.block.cast<T>(),     idx<4>({0,1,2,3},{0,2,3,1}));
    return result(0).real() / superblock->chain_length / superblock->chain_length;
}


double class_measurement::get_second_cumulant(){
    return first_moment();
}



double class_measurement::get_energy(){
    switch (sim){
        case SimulationType::iDMRG: case SimulationType::fDMRG: case SimulationType::FES_iDMRG:
            return get_expectationvalue(superblock->H.M);


        case SimulationType::iTEBD: case SimulationType::FES_iTEBD:
        {
            auto theta  = superblock->MPS.get_theta();

            Tensor<double,0> E1     = superblock->H.H_asTensor
                    .contract(theta.conjugate(), idx<2>({0,1},{0,1}))
                    .contract(theta,             idx<4>({0,1,2,3},{0,1,2,3}));
            return E1(0);
        }

        case SimulationType::NONE:
            std::cerr << "Simulation Type has not been set." << std::endl;
            exit(1);
    }
    return 0;

}

double class_measurement::get_entropy(){
    Tensor<Scalar,0> result1  = -superblock->MPS.LA.square()
            .contract(superblock->MPS.LA.square().log().eval(), idx<1>({0},{0}));
    return static_cast<double>(result1(0)) ;
}


double class_measurement::get_variance(){
    switch (sim) {
        case SimulationType::iDMRG: case SimulationType::fDMRG:  case SimulationType::FES_iDMRG:
            if (superblock->MPS.LA.size() == superblock->MPS.LB.size()) {
                return get_second_cumulant();
            }else{
                return 1;
            }


        case SimulationType::iTEBD: case SimulationType::FES_iTEBD:
        {
            Tensor<double, 4> theta = superblock->MPS.get_theta();
            Tensor<double, 4> H_sq = superblock->H.H_asTensor.contract(superblock->H.H_asTensor,
                                                                       idx<2>({0, 1}, {2, 3}));
            Tensor<double, 0> Var = H_sq.contract(theta.conjugate(), idx<2>({0, 1}, {0, 1}))
                    .contract(theta, idx<4>({0, 1, 2, 3}, {0, 1, 2, 3}));
            return (Var(0) - std::pow(get_energy(), 2));
        }

        case SimulationType::NONE:
            std::cerr << "Simulation Type has not been set." << std::endl;
            exit(1);
    }
    return 0;
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

