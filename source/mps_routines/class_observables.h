//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include <Eigen/Eigenvalues>
#include <complex>
#include <IO/class_custom_cout.h>
#include <general/class_eig_arpack_wrapper.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps.h>
#include <IO/class_multidata_buffer.h>
//using namespace std::complex_literals;

enum class SimulationType {iDMRG,fDMRG, FES_iDMRG, iTEBD, FES_iTEBD};
class class_observables {
public:
    using Scalar         = class_mps::Scalar;
private:
    class_superblock &superblock;
    class_multidata_buffer multibuffer;
    class_custom_cout ccout;
    SimulationType sim;

    using MapType = std::map<SimulationType, std::string>;
    MapType Sim2String;
    template <typename T>
    double find_slope(std::vector<T> &xvec, std::vector<T> &yvec ){
        Eigen::Array<T,Eigen::Dynamic, 1> X = Eigen::Map<Eigen::Array<T,Eigen::Dynamic, 1>>(xvec.data(), xvec.size());             //Cast to eigen array
        Eigen::Array<T,Eigen::Dynamic, 1> Y = Eigen::Map<Eigen::Array<T,Eigen::Dynamic, 1>>(yvec.data(), yvec.size());             //Cast to eigen array
        auto N = xvec.size();
        double slope = (N * X.cwiseProduct(Y).sum() - X.sum()*Y.sum() )/(N * X.cwiseProduct(X).sum() - X.sum()*X.sum() );
//        return ((X - X.mean()).cwiseProduct(Y - Y.mean())).sum() / fmax(1, (X - X.mean()).cwiseAbs2().sum());
        return slope;
    }

public:

    double first_moment(){

        superblock.update_bond_dimensions();
        if (superblock.MPS.LA.size() != superblock.MPS.LB.size()) { return 1.0;}
        double h = 0.0000001;

        double a = 0.0;
        double a1  = a - h;

        using T = std::complex<double>;
        class_eig_arpack eig;

        Textra::array<4> shape4 = superblock.shape4;
        long sizeL     = shape4[1] * shape4[1];
        long sizeR     = shape4[3] * shape4[3];

        Textra::Tensor<T,4> theta = superblock.MPS.thetaR().cast<T>();
        Textra::Tensor<T,3> A     = superblock.MPS.A().cast<T>();
        Textra::Tensor<T,2> transf_Ga1 = theta.contract(superblock.H.compute_G(a1), Textra::idx<2>({0,1},{0,1})).contract(theta.conjugate(), Textra::idx<2>({2,3},{0,1})).shuffle(Textra::array4{0,2,1,3}).reshape(Textra::array2{sizeL,sizeR}) ;


        Textra::Tensor<T,2> transf_ID = A.contract(A.conjugate(), Textra::idx<1>({0},{0}))
                                     .shuffle(Textra::array4{0,2,1,3})
                                     .reshape(Textra::array2{sizeL,sizeR}) ;

        Textra::Tensor<T,1> lambda_Ga1     = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Ga1, 1);
        T ga1  = (lambda_Ga1(0));


        double H2 = get_expectationvalue_sq(superblock.H.M);
        double E2 = get_expectationvalue(superblock.H.M) * get_expectationvalue(superblock.H.M);


        variance1 = H2-E2;
        variance2 = std::fabs((std::log(std::pow(std::fabs(ga1),2))));
        return std::log(std::pow(std::abs(ga1),2));
    }


    class_observables(class_superblock &superblockRef, SimulationType sim_);
    double variance1,variance2,variance3;
    double get_expectationvalue(const Textra::Tensor<double,4> &MPO);
    double get_expectationvalue(const Textra::Tensor<std::complex<double>,4> &MPO);
    double get_expectationvalue_sq(const Textra::Tensor<double,4> &MPO);
    double get_expectationvalue_sq(const Textra::Tensor<std::complex<double>,4> &MPO);
    double get_energy();                 /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
    double get_entropy();                /*! Computes the current entropy \f$ S = - \sum_n \lambda_n log( \lambda_n) \f$, where \f$\lambda_n \f$ are elements of \f$ \Lambda^A\f$ */
    double get_variance();               /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance1();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance2();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance3();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_truncation_error();
    double get_second_cumulant();
    long   get_chi_max();
    long   get_chi();
    long   get_chain_length();

    void print_status_full();                /*!< Print out status of all observables.*/
    void print_status_update(int step = 0);  /*!< Print out status of all observables.*/


};


#endif //DMRG_CLASS_OBSERVABLES_H
