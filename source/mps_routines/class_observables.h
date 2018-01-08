//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include <general/class_eig_arpack_wrapper.h>
#include <Eigen/Eigenvalues>
#include <complex>
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
    template <typename T>
    double find_slope(std::vector<T> &xvec, std::vector<T> &yvec ){
        Eigen::Array<T,Eigen::Dynamic, 1> X = Map<Eigen::Array<T,Eigen::Dynamic, 1>>(xvec.data(), xvec.size());             //Cast to eigen array
        Eigen::Array<T,Eigen::Dynamic, 1> Y = Map<Eigen::Array<T,Eigen::Dynamic, 1>>(yvec.data(), yvec.size());             //Cast to eigen array
        auto N = xvec.size();
        double slope = (N * X.cwiseProduct(Y).sum() - X.sum()*Y.sum() )/(N * X.cwiseProduct(X).sum() - X.sum()*X.sum() );
//        return ((X - X.mean()).cwiseProduct(Y - Y.mean())).sum() / fmax(1, (X - X.mean()).cwiseAbs2().sum());
        return slope;
    }

public:

    double first_moment(){

        superblock.update_bond_dimensions();
        if (superblock.MPS.LA.size() != superblock.MPS.LB.size()) { return 1.0;}
        double ha = 0.001;
        double hb = 0.00001;
        double hc = 0.0000001;

        double a = 0.0;
        double a0  = a;
        double a1  = a - ha;
        double a2  = a + ha;

        double b1  = a - hb;
        double b2  = a + hb;

        double c1  = a - hc;
        double c2  = a + hc;

        using T = std::complex<double>;


        Textra::array<4> shape4 = superblock.shape4;
        Textra::array<3> shape_Lblock = {shape4[1], shape4[1], superblock.H.M.dimension(0)};
        Textra::array<3> shape_Rblock = {shape4[3], shape4[3], superblock.H.M.dimension(1)};
        long sizeL     = shape4[1] * shape4[1];
        long sizeR     = shape4[3] * shape4[3];




        Tensor<T,4> theta = superblock.MPS.thetaR().cast<T>();
        Tensor<T,3> A     = superblock.MPS.A().cast<T>();
        Tensor<T,2> transf_G0 = theta.contract(superblock.H.compute_G(a0), idx<2>({0,1},{0,1})).contract(theta.conjugate() , idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Ga1 = theta.contract(superblock.H.compute_G(a1), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Ga2 = theta.contract(superblock.H.compute_G(a2), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Gb1 = theta.contract(superblock.H.compute_G(b1), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Gb2 = theta.contract(superblock.H.compute_G(b2), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Gc1 = theta.contract(superblock.H.compute_G(c1), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Gc2 = theta.contract(superblock.H.compute_G(c2), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;

        Tensor<T,2> transf_F0 = theta.contract(superblock.H.compute_F(a0), idx<2>({0,1},{0,1})).contract(theta.conjugate(), idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Fa1 = theta.contract(superblock.H.compute_F(a1), idx<2>({0,1},{0,1})).contract(theta.conjugate(),idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Fa2 = theta.contract(superblock.H.compute_F(a2), idx<2>({0,1},{0,1})).contract(theta.conjugate(),idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Fb1 = theta.contract(superblock.H.compute_F(b1), idx<2>({0,1},{0,1})).contract(theta.conjugate(),idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Fb2 = theta.contract(superblock.H.compute_F(b2), idx<2>({0,1},{0,1})).contract(theta.conjugate(),idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Fc1 = theta.contract(superblock.H.compute_F(c1), idx<2>({0,1},{0,1})).contract(theta.conjugate(),idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;
        Tensor<T,2> transf_Fc2 = theta.contract(superblock.H.compute_F(c2), idx<2>({0,1},{0,1})).contract(theta.conjugate(),idx<2>({2,3},{0,1})).shuffle(array4{0,2,1,3}).reshape(array2{sizeL,sizeR}) ;



        Tensor<T,2> transf_ID = A.contract(A.conjugate(), idx<1>({0},{0}))
                                     .shuffle(array4{0,2,1,3})
                                     .reshape(array2{sizeL,sizeR}) ;


        class_eig_arpack eig;
        Tensor<T,1> lambda_G0      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_G0, 1);
        Tensor<T,1> lambda_Ga1      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Ga1, 1);
        Tensor<T,1> lambda_Ga2      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Ga2, 1);
        Tensor<T,1> lambda_Gb1      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Gb1, 1);
        Tensor<T,1> lambda_Gb2      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Gb2, 1);
        Tensor<T,1> lambda_Gc1      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Gc1, 1);
        Tensor<T,1> lambda_Gc2      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Gc2, 1);
        Tensor<T,1> lambda_F0      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_F0, 1);
        Tensor<T,1> lambda_Fa1      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Fa1, 1);
        Tensor<T,1> lambda_Fa2      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Fa2, 1);
        Tensor<T,1> lambda_Fb1      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Fb1, 1);
        Tensor<T,1> lambda_Fb2      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Fb2, 1);
        Tensor<T,1> lambda_Fc1      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Fc1, 1);
        Tensor<T,1> lambda_Fc2      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_Fc2, 1);
        Tensor<T,1> lambda_ID_R    = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transf_ID, 1);
        Tensor<T,1> lambda_ID_L    = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::L, false>(transf_ID, 1);

//        T g0  = std::pow(lambda_G0(0)/lambda_ID_R(0), 0.5);
//        T g1  = std::pow(lambda_G1(0)/lambda_ID_R(0), 0.5);
//        T g2  = std::pow(lambda_G2(0)/lambda_ID_R(0), 0.5);
//        T f0  = std::pow(lambda_F0(0)/lambda_ID_R(0), 0.5);
//        T f1  = std::pow(Lambda_F1(0)/lambda_ID_R(0), 0.5);
//        T f2  = std::pow(lambda_F2(0)/lambda_ID_R(0), 0.5);


        cout << "lambda_G0    " <<lambda_G0   << endl;
        cout << "lambda_Ga1   " <<lambda_Ga1  << endl;
        cout << "lambda_Ga2   " <<lambda_Ga2  << endl;
        cout << "lambda_Gb1   " <<lambda_Gb1  << endl;
        cout << "lambda_Gb2   " <<lambda_Gb2  << endl;
        cout << "lambda_Gc1   " <<lambda_Gc1  << endl;
        cout << "lambda_Gc2   " <<lambda_Gc2  << endl;

        cout << "lambda_F0    " <<lambda_F0   << endl;
        cout << "lambda_Fa1   " <<lambda_Fa1  << endl;
        cout << "lambda_Fa2   " <<lambda_Fa2  << endl;
        cout << "lambda_Fb1   " <<lambda_Fb1  << endl;
        cout << "lambda_Fb2   " <<lambda_Fb2  << endl;
        cout << "lambda_Fc1   " <<lambda_Fc1  << endl;
        cout << "lambda_Fc2   " <<lambda_Fc2  << endl;
        cout << "lambda_ID_R  " <<lambda_ID_R << endl;
        cout << "lambda_ID_L  " <<lambda_ID_L << endl;

        T g0   = lambda_G0(0);
        T ga1  = (lambda_Ga1(0));
        T ga2  = (lambda_Ga2(0));
        T gb1  = (lambda_Gb1(0));
        T gb2  = (lambda_Gb2(0));
        T gc1  = (lambda_Gc1(0));
        T gc2  = (lambda_Gc2(0));

        T f0   = lambda_F0(0);
        T fa1  = (lambda_Fa1(0));
        T fa2  = (lambda_Fa2(0));
        T fb1  = (lambda_Fb1(0));
        T fb2  = (lambda_Fb2(0));
        T fc1  = (lambda_Fc1(0));
        T fc2  = (lambda_Fc2(0));



        T gpa  =  (ga2-ga1)/2.0/ha/std::complex<double>(0,1);
        T gpb  =  (gb2-gb1)/2.0/hb/std::complex<double>(0,1);
        T gpc  =  (gc2-gc1)/2.0/hc/std::complex<double>(0,1);


        T gppa = (ga2 - 2.0*g0 + ga1)/(ha*ha);
        T gppb = (gb2 - 2.0*g0 + gb1)/(hb*hb);
        T gppc = (gc2 - 2.0*g0 + gc1)/(hc*hc);




        T fpa  = (fa2-fa1)/2.0/ha;
        T fpb  = (fb2-fb1)/2.0/hb;
        T fpc  = (fc2-fc1)/2.0/hc;
        T fppa = (fa2 - 2.0 + fa1)/(ha*ha);
        T fppb = (fb2 - 2.0 + fb1)/(hb*hb);
        T fppc = (fc2 - 2.0 + fc1)/(hc*hc);


        T lfpa = (log(fa2)-log(fa1))/2.0/ha;
        T lfpb = (log(fb2)-log(fb1))/2.0/hb;
        T lfpc = (log(fc2)-log(fc1))/2.0/hc;


        double H2 = get_expectationvalue_sq(superblock.H.M);
        double E2 = get_expectationvalue(superblock.H.M) * get_expectationvalue(superblock.H.M);

        cout << setprecision(16);
        cout << endl;

        double E = -1.2732395447351625;

        cout << "G0          : "   << g0 << endl;
        cout << "Ga1         : "   << ga1 << endl;
        cout << "Gb1         : "   << gb1 << endl;
        cout << "Gc1         : "   << gc1 << endl;
        cout << "Ga2         : "   << ga2 << endl;
        cout << "Gb2         : "   << gb2 << endl;
        cout << "Gc2         : "   << gc2 << endl;

        cout << "F0          : "   << f0 << endl;
        cout << "Fa1         : "   << fa1 << endl;
        cout << "Fb1         : "   << fb1 << endl;
        cout << "Fc1         : "   << fc1 << endl;
        cout << "Fa2         : "   << fa2 << endl;
        cout << "Fb2         : "   << fb2 << endl;
        cout << "Fc2         : "   << fc2 << endl;

        cout << "Energy     F'(a)                    : " <<right << fpa << " Ea - E: " << fpa - E <<  endl;
        cout << "Energy     F'(b)                    : " <<right << fpb << " Eb - E: " << fpb - E <<  endl;
        cout << "Energy     F'(c)                    : " <<right << fpc << " Ec - E: " << fpc - E <<  endl;
        cout << "Variance   <M2>                     : " <<right << (fa2.real() + fa1.real() - 2)/(ha*ha)  << endl;
        cout << "Variance   <M2> - E²                : " <<right << fppa   << endl;
        cout << "Variance   F''(a) - E²              : " <<right << fppa   << endl;
        cout << "Variance   F''(b) - E²              : " <<right << fppb   << endl;
        cout << "Variance   F''(c) - E²              : " <<right << fppc   << endl;
        cout << "Variance   log F(a) + log F(-a))    : " <<right <<  std::log(fa2) + std::log(fa1)  << endl;
        cout << "Variance   log F(b) + log F(-b))    : " <<right <<  std::log(fb2) + std::log(fb1)  << endl;
        cout << "Variance   log F(c) + log F(-c))    : " <<right <<  std::log(fc2) + std::log(fc1)  << endl;
        cout << '\n';

        cout << "Energy     G'(a)                    : " <<right << gpa << " Ea - E: " << gpa - E <<  endl;
        cout << "Energy     G'(b)                    : " <<right << gpb << " Eb - E: " << gpb - E <<  endl;
        cout << "Energy     G'(c)                    : " <<right << gpc << " Ec - E: " << gpc - E <<  endl;
        cout << "Variance   G''(a) - E²              : " <<right << gppa << " - " << std::norm(gpa) << " = " << gppa - std::norm(gpa)  << endl;
        cout << "Variance   G''(b) - E²              : " <<right << gppb << " - " << std::norm(gpb) << " = " << gppb - std::norm(gpb)  << endl;
        cout << "Variance   G''(c) - E²              : " <<right << gppc << " - " << std::norm(gpc) << " = " << gppc - std::norm(gpc)  << endl;

        cout << "Variance   log G(a) + log G(-a)     : " <<right << std::log(ga2) + std::log(ga1) << endl;
        cout << "Variance   log G(b) + log G(-b)     : " <<right << std::log(gb2) + std::log(gb1) << endl;
        cout << "Variance   log G(c) + log G(-c)     : " <<right << std::log(gc2) + std::log(gc1) << endl;

        cout << "Variance   log(|G(a)|²)     >       : " <<right << std::log(std::pow(std::abs(ga2),2)) << endl;
        cout << "Variance   log(|G(b)|²)     >       : " <<right << std::log(std::pow(std::abs(gb2),2)) << endl;
        cout << "Variance   log(|G(c)|²)     >       : " <<right << std::log(std::pow(std::abs(gc2),2)) << endl;




        cout << '\n';
        cout << "Variance   2 MPO's                  : " <<  H2 - E2 << " H2: " << H2 << " E2: " << E2 <<  endl;
//        if (H2-E2 < 0) {exit(1);}

        variance1 = H2-E2;
        variance2 = std::fabs((std::log(std::pow(std::fabs(gb2),2))));
        variance3 = std::fabs((std::log(std::pow(std::fabs(gc2),2))));
        return std::log(std::pow(std::abs(gc2),2));
    }


    class_observables(class_superblock &superblockRef, SimulationType sim_);
    double variance1,variance2,variance3;
    double get_expectationvalue(const Tensor<double,4> &MPO);
    double get_expectationvalue(const Tensor<std::complex<double>,4> &MPO);
    double get_expectationvalue_sq(const Tensor<double,4> &MPO);
    double get_expectationvalue_sq(const Tensor<std::complex<double>,4> &MPO);
    double get_energy();                /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
    double get_entropy();               /*! Computes the current entropy \f$ S = - \sum_n \lambda_n log( \lambda_n) \f$, where \f$\lambda_n \f$ are elements of \f$ \Lambda^A\f$ */
    double get_variance();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance1();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance2();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance3();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_truncation_error();
    double get_second_cumulant();
    long   get_chi_max();
    long   get_chi();
    long   get_chain_length();

    void print_status_full(int verbosity);   /*!< Print out status of all observables.*/
    void print_status_update(int step = 0);  /*!< Print out status of all observables.*/


};


#endif //DMRG_CLASS_OBSERVABLES_H
