//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include <general/class_eig_arpack_wrapper.h>
#include <Eigen/Eigenvalues>
#include "class_superblock.h"
using namespace std::complex_literals;

enum class SimulationType {iDMRG,fDMRG, FES_iDMRG, iTEBD, FES_iTEBD};
class class_observables {
public:
    using Scalar         = class_mps::Scalar;
private:
    class_superblock &superblock;
    SimulationType sim;
    using MapType = std::map<SimulationType, std::string>;
    MapType Sim2String;
//    Tensor<std::complex<double>,4> get_Lblock_sq ();
//    Tensor<std::complex<double>,4> get_Rblock_sq ();
public:

    double first_moment(){

        superblock.update_bond_dimensions();
        if (superblock.MPS.LA.size() != superblock.MPS.LB.size()) { return 1.0;}
        double h = 0.000001;
        double a = 0.0;
        double a0 = a;
        double a1 = a - h;
        double a2 = a + h;
        using T = std::complex<double>;


        Textra::array<4> shape4 = superblock.shape4;
        Textra::array<3> shape_Lblock = {shape4[1], shape4[1], superblock.H.M.dimension(0)};
        Textra::array<3> shape_Rblock = {shape4[3], shape4[3], superblock.H.M.dimension(1)};
        long sizeL     = shape4[1] * shape4[1];
        long sizeR     = shape4[3] * shape4[3];




        Tensor<T,4> theta = superblock.MPS.thetaR().cast<T>();
        Tensor<T,2> transf_G = theta.contract(superblock.H.compute_G(h), idx<2>({0,1},{0,1}))
                .contract(theta.conjugate(),          idx<2>({2,3},{0,1}))
                .shuffle(array4{0,2,1,3})
                .reshape(array2{sizeL,sizeR}) ;

        Tensor<T,2> transf_F0 = theta.contract(superblock.H.compute_F(a0), idx<2>({0,1},{0,1}))
                .contract(theta.conjugate(),          idx<2>({2,3},{0,1}))
                .shuffle(array4{0,2,1,3})
                .reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_F1 = theta.contract(superblock.H.compute_F(a1), idx<2>({0,1},{0,1}))
                                       .contract(theta.conjugate(),          idx<2>({2,3},{0,1}))
                                       .shuffle(array4{0,2,1,3})
                                       .reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_F2 = theta.contract(superblock.H.compute_F(a2), idx<2>({0,1},{0,1}))
                                     .contract(theta.conjugate(),          idx<2>({2,3},{0,1}))
                                     .shuffle(array4{0,2,1,3})
                                     .reshape(array2{sizeL,sizeR}) ;

        Tensor<T,2> transf_ID = theta.contract(theta.conjugate(), idx<2>({0,1},{0,1}))
                                          .shuffle(array4{0,2,1,3})
                                          .reshape(array2{sizeL,sizeR}) ;


        class_eig_arpack eig;
        auto lambda_G = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_G , 1);
        auto lambda_0 = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_F0, 1);
        auto lambda_1 = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_F1, 1);
        auto lambda_2 = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_F2, 1);
        auto lambda_R = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_ID, 1);
        auto lambda_L = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::L, false>(transf_ID, 1);
        cout << setprecision(16);
        cout << endl;
        cout << "lambda_G: " << lambda_G(0) << endl;
        cout << "lambda_0: " << lambda_0(0) << endl;
        cout << "lambda_1: " << lambda_1(0) << endl;
        cout << "lambda_2: " << lambda_2(0) << endl;
        cout << "lambda_R: " << lambda_R(0) << endl;
        cout << "lambda_L: " << lambda_L(0) << endl;
        double g  = std::pow(std::norm(lambda_G(0))/std::norm(lambda_R(0)), 0.5);
        double f0 = std::pow(std::norm(lambda_0(0))/std::norm(lambda_R(0)), 0.5);
        double f1 = std::pow(std::norm(lambda_1(0))/std::norm(lambda_R(0)), 0.5);
        double f2 = std::pow(std::norm(lambda_2(0))/std::norm(lambda_R(0)), 0.5);
        double energy   = (f2-f1)/2.0/h;
        double variance = (f2 - f0*2.0 + f1)/(h*h) - energy*energy;
        cout << "Energy: " << energy << endl;
        cout << "Variance: " << variance  << endl;
        cout << "Variance: " << (std::log(f2) + std::log(f1))/(h*h)  << endl;
        cout << "Variance: " << std::log(std::pow(g,2)) << endl;
        return (f2-f1)/2.0/h;
//        return  std::pow((lambda_2(0) - lambda_1(0))/lambda_R(0) , 0.5)/(2*h);
//        return  (std::norm(lambda_2(0)) - std::norm(lambda_1(0)))/(2*h)/std::norm(lambda_R(0));

    }


    class_observables(class_superblock &superblockRef, SimulationType sim_);

    double get_expectationvalue(const Tensor<double,4> &MPO);
    double get_expectationvalue(const Tensor<std::complex<double>,4> &MPO);

    double get_energy();                /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
    double get_entropy();               /*! Computes the current entropy \f$ S = - \sum_n \lambda_n log( \lambda_n) \f$, where \f$\lambda_n \f$ are elements of \f$ \Lambda^A\f$ */
    double get_variance();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_truncation_error();
    double get_second_cumulant();
    long   get_chi_max();
    long   get_chi();
    long   get_chain_length();

    void print_status_full(int verbosity);   /*!< Print out status of all observables.*/
    void print_status_update(int step = 0);  /*!< Print out status of all observables.*/


};


#endif //DMRG_CLASS_OBSERVABLES_H
