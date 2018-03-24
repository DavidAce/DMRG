//
// Created by david on 2018-02-21.
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
#include <general/class_svd_wrapper.h>
#include <Eigen/QR>

using namespace std;
using namespace Textra;
using Scalar = class_measurement::Scalar;
using namespace std::complex_literals;


Scalar class_measurement::moment_generating_function(std::shared_ptr<class_mps> MPS_original,
                                                               std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
    std::shared_ptr<class_mps> MPS_evolved = std::make_shared<class_mps>(*MPS_original);
    class_SVD<Scalar> SVD;
    class_arpackpp_wrapper eig;
    SVD.setThreshold(1e-12);
    eig.setThreshold(1e-12);
    for (auto &Op: Op_vec) {
        //Evolve
        Textra::Tensor<Scalar, 4> theta_evo = Op.contract(MPS_evolved->get_theta(), idx({0, 1}, {0, 2})).shuffle(array4{0, 2, 1, 3});
        auto chiL = theta_evo.dimension(1);
        auto chiR = theta_evo.dimension(3);
        auto[U, S, V] = SVD.schmidt(theta_evo, 2, chiL, chiR);
        MPS_evolved->LA= S;
        MPS_evolved->GA = asDiagonalInversed(MPS_evolved->LB).contract(U, idx({1}, {1})).shuffle(array3{1, 0, 2});
        MPS_evolved->GB = V.contract(asDiagonalInversed(MPS_evolved->LB), idx({2}, {0}));
        if (&Op != &Op_vec.back()) {
            MPS_evolved->swap_AB();
        }

    }
    long sizeLB = MPS_evolved->LB.size() * MPS_evolved->LB.size();

    //Normalize
    Tensor<Scalar,2> theta_evn_mat       = MPS_evolved->get_transfer_matrix_theta_evn().reshape(array2{sizeLB,sizeLB});
    VectorType<Scalar> eigval_R_evn      = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(theta_evn_mat.data(),sizeLB,sizeLB, 1);
    auto new_theta_evn_normalized        = MPS_evolved->get_theta_evn(sqrt(eigval_R_evn(0)));

    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original->theta_evn_normalized.dimension(1);
    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original->theta_evn_normalized.dimension(3);

    Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
            .contract(MPS_original->theta_evn_normalized.conjugate(), idx({0,2},{0,2}))
            .shuffle(array4{0,2,1,3})
            .reshape(array2{sizeL,sizeR});
    //Compute the characteristic function G(a).
    eig.setThreshold(1e-14);
    VectorType <Scalar> lambdaG  = eig.solve_dominant<arpack::Form::COMPLEX, arpack::Ritz::LM, arpack::Side::R, false>(transfer_matrix_G.data(), transfer_matrix_G.dimension(0),transfer_matrix_G.dimension(1), 1);
    return lambdaG(0);

}



std::pair<double,double> class_measurement::compute_infinite_moments_G(Scalar a, std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
    using T = Scalar;
    using LT = std::complex<long double>;
    LT al = (LT) a;
    //The following only works if superblock->MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    LT lambdaG  = (LT) moment_generating_function(superblock->MPS, Op_vec);
    LT l        = (LT) (superblock->H->mps_sites*1.0l);
    LT G        = pow(lambdaG,1.0l/l);
    LT logG     = log(lambdaG) * (1.0l/l);
    LT logGc    = log(conj(lambdaG) ) * (1.0l/l);
    LT O        = (logG - logGc)/(2.0l*al);
    LT VarO     = 2.0l*log(abs(G))/ (al*al);
    return std::make_pair(real(O), real(VarO));




//
//    T l2 = 2.0+0.0i;
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esA(Tensor2_to_Matrix(transf_GA));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esB(Tensor2_to_Matrix(transf_GB));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esC(Tensor2_to_Matrix(transf_GC));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esD(Tensor2_to_Matrix(transf_GD));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esiA(Tensor2_to_Matrix(transf_ID_A));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esiB(Tensor2_to_Matrix(transf_ID_B));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esiC(Tensor2_to_Matrix(transf_ID_C));
//    Eigen::ComplexEigenSolver<Textra::MatrixType<T>> esiD(Tensor2_to_Matrix(transf_ID_D));
//    Eigen::ArrayXcd eigvA  = esA.eigenvalues();
//    Eigen::ArrayXcd eigvB  = esB.eigenvalues();
//    Eigen::ArrayXcd eigvC  = esC.eigenvalues();
//    Eigen::ArrayXcd eigvD  = esD.eigenvalues();
//    Eigen::ArrayXcd eigviA = esiA.eigenvalues();
//    Eigen::ArrayXcd eigviB = esiB.eigenvalues();
//    Eigen::ArrayXcd eigviC = esiC.eigenvalues();
//    Eigen::ArrayXcd eigviD = esiD.eigenvalues();
////    Eigen::ArrayXcd eigviB = esib.eigenvalues();
//    T largestA = eigvA.tail(1).coeff(0);
//    T largestB = eigvB.tail(1).coeff(0);
//    T largestC = eigvC.tail(1).coeff(0);
//    T largestD = eigvD.tail(1).coeff(0);
//    T largestia = eigviA.tail(1).coeff(0);
//    T largestib = eigviB.tail(1).coeff(0);
//    T largestic = eigviC.tail(1).coeff(0);
//    T largestid = eigviD.tail(1).coeff(0);
////    T largestib = eigviB.tail(1).coeff(0);
//    Eigen::ArrayXcd eigvA_j = eigvA.head(eigvA.size()-1);
//    Eigen::ArrayXcd eigvB_j = eigvB.head(eigvB.size()-1);
//    Eigen::ArrayXcd eigvC_j = eigvC.head(eigvC.size()-1);
//    Eigen::ArrayXcd eigvD_j = eigvD.head(eigvD.size()-1);
//    Eigen::ArrayXcd eigviA_j = eigviA.head(eigviA.size()-1);
//    Eigen::ArrayXcd eigviB_j = eigviB.head(eigviB.size()-1);
//    Eigen::ArrayXcd eigviC_j = eigviC.head(eigviC.size()-1);
//    Eigen::ArrayXcd eigviD_j = eigviD.head(eigviD.size()-1);
////    Eigen::ArrayXcd eigviB_j = eigviB.head(eigviB.size()-1);
//    T logg2a  = (L/l2) * std::log(largestA) + std::log(1.0 + (L/l2* (eigvA_j/largestA).log()).exp().sum());
//    T logg2b  = (L/l2) * std::log(largestB) + std::log(1.0 + (L/l2* (eigvB_j/largestB).log()).exp().sum());
//    T logg2c  = (L/l2) * std::log(largestC) + std::log(1.0 + (L/l2* (eigvC_j/largestC).log()).exp().sum());
//    T logg2d  = (L/l2) * std::log(largestD) + std::log(1.0 + (L/l2* (eigvD_j/largestD).log()).exp().sum());
//    T logg2ac = (L/l2) * std::log(conj(largestA)) + std::log(1.0 + (L/l2* (eigvB_j/largestA).conjugate().log()).exp().sum());
//    T logg2bc = (L/l2) * std::log(conj(largestB)) + std::log(1.0 + (L/l2* (eigvB_j/largestB).conjugate().log()).exp().sum());
//    T logg2cc = (L/l2) * std::log(conj(largestC)) + std::log(1.0 + (L/l2* (eigvC_j/largestC).conjugate().log()).exp().sum());
//    T logg2cd = (L/l2) * std::log(conj(largestD)) + std::log(1.0 + (L/l2* (eigvD_j/largestD).conjugate().log()).exp().sum());
//    T logg2ia = (L/l2) * std::log(largestia) + std::log(1.0 + (L/l2* (eigviA_j/largestia).log()).exp().sum());
//    T logg2ib = (L/l2) * std::log(largestib) + std::log(1.0 + (L/l2* (eigviB_j/largestib).log()).exp().sum());
//    T logg2ic = (L/l2) * std::log(largestic) + std::log(1.0 + (L/l2* (eigviC_j/largestic).log()).exp().sum());
//    T logg2id = (L/l2) * std::log(largestid) + std::log(1.0 + (L/l2* (eigviD_j/largestid).log()).exp().sum());
//    T logG2A = logg2a  - logg2ia;
//    T logG2B = logg2b  - logg2ib;
//    T logG2C = logg2c  - logg2ic;
//    T logG2D = logg2d  - logg2id;
//    T logG2Ac= logg2ac - logg2ia;
//    T logG2Bc= logg2bc - logg2ib;
//    T logG2Cc= logg2cc - logg2ic;
//    T logG2Dc= logg2cd - logg2id;
//    T g2a  = eigvA.pow(L/l2).sum();
//    T g2b  = eigvB.pow(L/l2).sum();
//    T g2c  = eigvC.pow(L/l2).sum();
//    T g2d  = eigvD.pow(L/l2).sum();
//    T g2ia = eigviA.pow(L/l2).sum();
//    T g2ib = eigviB.pow(L/l2).sum();
//    T g2ic = eigviC.pow(L/l2).sum();
//    T g2id = eigviD.pow(L/l2).sum();
//    T G2A     = g2a/g2ia;
//    T G2B     = g2b/g2ib;
//    T G2C     = g2c/g2ic;
//    T G2D     = g2d/g2id;

//    T e2a     = (logG2A - logG2Ac)/(2.0*a);
//    T e2b     = (logG2B - logG2Bc)/(2.0*b);
//    T e2c     = (logG2C - logG2Cc)/(2.0*c);
//    T e2d     = (logG2D - logG2Dc)/(2.0*d);
//    T varHE2a  = (logG2A + logG2Ac)/(a*a)/L;
//    T varHE2b  = (logG2B + logG2Bc)/(b*b)/L;
//    T varHE2c  = (logG2C + logG2Cc)/(c*c)/L;
//    T varHE2d  = (logG2D + logG2Dc)/(d*d)/L;
//    T test_varHE2a  =  (logG2A + logG2Ac)/ (a*a*L)  ;
//    T test_varHE2b  =  (logG2B + logG2Bc)/ (b*b*L)  ;
//    T test_varHE2c  =  (logG2C + logG2Cc)/ (c*c*L)  ;
//    T test_varHE2d  =  (logG2D + logG2Dc)/ (d*d*L)  ;
//    std::cout << "E2A           = " << setw(30)  << e2a/L                                               << std::endl;
//    std::cout << "E2B           = " << setw(30)  << e2b/L                                               << std::endl;
//    std::cout << "E2C           = " << setw(30)  << e2c/L                                               << std::endl;
//    std::cout << "E2D           = " << setw(30)  << e2d/L                                               << std::endl;
//    std::cout << "varHE2a       = " << setw(45)  << varHE2a           << " -> " <<  std::log10(varHE2a )   << std::endl;
//    std::cout << "varHE2b       = " << setw(45)  << varHE2b           << " -> " <<  std::log10(varHE2b )   << std::endl;
//    std::cout << "varHE2c       = " << setw(45)  << varHE2c           << " -> " <<  std::log10(varHE2c )   << std::endl;
//    std::cout << "varHE2d       = " << setw(45)  << varHE2d           << " -> " <<  std::log10(varHE2d )   << std::endl;
//    std::cout << "test_varHE2a  = " << setw(45)  << test_varHE2a      << " -> " <<  std::log10(test_varHE2a)   << std::endl;
//    std::cout << "test_varHE2b  = " << setw(45)  << test_varHE2b      << " -> " <<  std::log10(test_varHE2b)   << std::endl;
//    std::cout << "test_varHE2c  = " << setw(45)  << test_varHE2c      << " -> " <<  std::log10(test_varHE2c)   << std::endl;
//    std::cout << "test_varHE2d  = " << setw(45)  << test_varHE2d      << " -> " <<  std::log10(test_varHE2d)   << std::endl;









    //    T G2A = g2a/std::abs(g2a);
//    T G2B = g2b/std::abs(g2b);
//    G2A /= abs(G2A);
//    G2B /= abs(G2B);
//    T e10     = (G1B - std::conj(G1B))/(2.0*a) * l ;
//    T e11     = (std::log(G1B) - std::log(std::conj(G1B)))/(2.0*a)*l;
//    T e1b     = (logG1B - conj(logG1B))/(2.0*a)*l;

//    T H10     = (G1B + std::conj(G1B) - 2.0) / (a*a)*l*l;
//    T H11     = (G1B + std::conj(G1B) - 2.0) / (a*a)*l*l;
//    T E10     = e10*std::conj(e10);
//    T E11     = e11*std::conj(e11);
//    T E12     = e1b*std::conj(e1b);

//    T varHE10 = H10 - E10;
//    T varHE11 = H10 - E11;
//    T varHE12 = std::log(std::pow(G1B*std::conj(G1B), 1.0/l))/(a*a);


//    T e20      = L * ( l_j_Ll1.cwiseProduct(l_j_p)  ).sum() / l_j_Ll.sum();
//    T e21     = (std::log(G2B) - std::log(conj(G2B)))/(2.0*a)*l2 ;
//    T e2b     = (logG2B - conj(logG2B))/(2.0*a)*l2;

//    T H20     = (G2B + std::conj(G2B) - 2.0) / (a*a*L*L)*l2*l2;
//    T E20     = (e20)* std::conj(e20) / L / L;
//    T E21     = (e21)*std::conj(e21) /L / L;
//    T E22     = (e2b)*std::conj(e2b) /L / L;
//    T varHE20 = (H20 - E20);
//    T varHE21 = (H20 - E21);
//    T varHE22 = (  std::log(G2B) + std::log( std::conj(G2B) )  )/(L*a*a);

}


double class_measurement::compute_infinite_variance_MPO(){
    Tensor<Scalar, 0> E2 =
            superblock->Lblock2->block
                    .contract(asDiagonal(superblock->MPS->LA), idx({0},{0}))
                    .contract(asDiagonal(superblock->MPS->LA), idx({0},{0}))
            .contract(superblock->Rblock2->block, idx({2, 3, 0, 1}, {0, 1, 2, 3}));
    double L = superblock->chain_length;
    return std::real( E2(0)/ (L - 2.0));
}

double class_measurement::compute_infinite_variance_H(){
    const Tensor<Scalar,4> & theta                      = superblock->MPS->theta;
    const Tensor<Scalar,4> & theta_sw                   = superblock->MPS->theta_sw;
    const Tensor<Scalar,4> & theta_evn                  = superblock->MPS->theta_evn_normalized;
    const Tensor<Scalar,4> & theta_odd                  = superblock->MPS->theta_odd_normalized;
    const Tensor<Scalar,3> & LBGA                       = superblock->MPS->LBGA;
    const Tensor<Scalar,3> & LAGB                       = superblock->MPS->LAGB;
    const Tensor<Scalar,2> & l_evn                      = superblock->MPS->l_evn;
    const Tensor<Scalar,2> & r_evn                      = superblock->MPS->r_evn;
    const Tensor<Scalar,2> & l_odd                      = superblock->MPS->l_odd;
    const Tensor<Scalar,2> & r_odd                      = superblock->MPS->r_odd;
    const Tensor<Scalar,4> & transfer_matrix_evn        = superblock->MPS-> transfer_matrix_evn;
    const Tensor<Scalar,4> & transfer_matrix_odd        = superblock->MPS-> transfer_matrix_odd;
    const Tensor<Scalar,4> & transfer_matrix_LBGA       = superblock->MPS-> transfer_matrix_LBGA;
    const Tensor<Scalar,4> & transfer_matrix_LAGB       = superblock->MPS-> transfer_matrix_LAGB;

//    Tensor<Scalar,4> h0_   = Textra::Matrix_to_Tensor(superblock->H->h[0], 2,2,2,2);
//    Tensor<Scalar,4> h1_   = Textra::Matrix_to_Tensor(superblock->H->h[1], 2,2,2,2);
//
//    Tensor<Scalar,0> EAB   =
//            theta_evn
//                    .contract(h0_,                   idx({0, 2}, {0, 1}))
//                    .contract(theta_evn.conjugate(),idx({2, 3}, {0, 2}))
//                    .contract(l_evn,                idx({0, 2}, {0, 1}))
//                    .contract(r_evn,                idx({0, 1}, {0, 1}));
//
//    Tensor<Scalar,0> EBA   =
//            theta_odd
//                    .contract(h1_, idx({0, 2}, {0, 1}))
//                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0,2}))
//                    .contract(l_odd, idx({0,2}, {0,1}))
//                    .contract(r_odd, idx({0,1},{0,1}));
////
//            Tensor<Scalar,4> I_EAB = Matrix_to_Tensor((EAB(0)*MatrixType<Scalar>::Identity(4,4)).eval() , 2,2,2,2);
//            Tensor<Scalar,4> I_EBA = Matrix_to_Tensor((EBA(0)*MatrixType<Scalar>::Identity(4,4)).eval() , 2,2,2,2);
//            Tensor<Scalar,4> h0   = h0_ - I_EAB;
//            Tensor<Scalar,4> h1   = h1_ - I_EBA;

    Tensor<Scalar,4> h0 =  Matrix_to_Tensor((superblock->H->h[0] - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
    Tensor<Scalar,4> h1 =  Matrix_to_Tensor((superblock->H->h[1] - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);

    Tensor<Scalar,0> E2AB =
             theta_evn
            .contract(h0,                     idx({0, 2}, {0, 1}))
            .contract(h0,                     idx({2, 3}, {0, 1}))
            .contract(theta_evn.conjugate(),  idx({2, 3}, {0, 2}))
            .contract(l_evn,                  idx({0, 2}, {0, 1}))
            .contract(r_evn,                  idx({0, 1}, {0, 1}));


    Tensor<Scalar, 0> E2BA =
             theta_odd
            .contract(h1,                    idx({0, 2}, {0, 1}))
            .contract(h1,                    idx({2, 3}, {0, 1}))
            .contract(theta_odd.conjugate(), idx({2, 3}, {0, 2}))
            .contract(l_odd,                 idx({0, 2}, {0, 1}))
            .contract(r_odd,                 idx({0, 1}, {0, 1}));



    Tensor<Scalar,5> thetaABA = theta_evn.contract(LBGA, idx({3},{1}));
    Tensor<Scalar,5> thetaBAB = theta_odd.contract(LAGB, idx({3},{1}));

    Tensor<Scalar,0> E2ABA_1  =
                      thetaABA
                     .contract(h1,                   idx({2,3},{0,1}))
                     .contract(h0,                   idx({0,3},{0,1}))
                     .contract(thetaABA.conjugate(), idx({3,4,2},{0,2,3}))
                     .contract(l_evn,                idx({0,2},{0,1}))
                     .contract(r_odd,                idx({0,1},{0,1})) ;

    Tensor<Scalar,0> E2BAB_1  =
                    thetaBAB
                    .contract(h1,                   idx({0,2},{0,1}))
                    .contract(h0,                   idx({4,1},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({2,3,4},{0,2,3}))
                    .contract(l_odd,                idx({0,2},{0,1}))
                    .contract(r_evn,                idx({0,1},{0,1})) ;

    Tensor<Scalar,0> E2ABA_2  =
                    thetaABA
                    .contract(h0,                   idx({0,2},{0,1}))
                    .contract(h1,                   idx({4,1},{0,1}))
                    .contract(thetaABA.conjugate(), idx({2,3,4},{0,2,3}))
                    .contract(l_evn,                idx({0,2},{0,1}))
                    .contract(r_odd,                idx({0,1},{0,1})) ;

    Tensor<Scalar,0> E2BAB_2  =
                     thetaBAB
                    .contract(h0,                   idx({2,3},{0,1}))
                    .contract(h1,                   idx({0,3},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({3,4,2},{0,2,3}))
                    .contract(l_odd,                idx({0,2},{0,1}))
                    .contract(r_evn,                idx({0,1},{0,1})) ;


    Tensor<Scalar,2> E2d_L_evn =
            theta_evn
                    .contract(h0,                    idx({0, 2}, {0, 1}))
                    .contract(theta_evn.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(l_evn,                 idx({0, 2}, {0, 1}));

    Tensor<Scalar,2> E2d_R_evn =
                     theta_evn
                    .contract(h0,                    idx({0, 2}, {0, 1}))
                    .contract(theta_evn.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(r_evn,                 idx({1, 3}, {0, 1}));

    Tensor<Scalar,2> E2d_L_odd  =
                     theta_odd
                    .contract(h1,                     idx({0, 2}, {0, 1}))
                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(l_odd,                  idx({0, 2}, {0, 1}));


    Tensor<Scalar,2> E2d_R_odd =
                     theta_odd
                    .contract(h1,                     idx({0, 2}, {0, 1}))
                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(r_odd,                  idx({1, 3}, {0, 1}));

    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
    Tensor<Scalar,4> fixpoint_evn = r_evn.contract(l_evn, pair);
    Tensor<Scalar,4> fixpoint_odd = r_odd.contract(l_odd, pair);

    long sizeLA = superblock->MPS->LA.size();
    long sizeLB = superblock->MPS->LB.size();
    Tensor<Scalar,2> one_minus_transfer_matrix_evn = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLB*sizeLB, sizeLA*sizeLA).eval()) - (transfer_matrix_evn-fixpoint_evn).reshape(array2{sizeLB*sizeLB, sizeLA*sizeLA});
    Tensor<Scalar,2> one_minus_transfer_matrix_odd = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLA*sizeLA, sizeLB*sizeLB).eval()) - (transfer_matrix_odd-fixpoint_odd).reshape(array2{sizeLA*sizeLA, sizeLB*sizeLB});

    class_SVD<Scalar> SVD;
    Tensor<Scalar,4> E_evn_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_evn).reshape(array4{sizeLB,sizeLB,sizeLA,sizeLA});
    Tensor<Scalar,4> E_odd_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_odd).reshape(array4{sizeLA,sizeLA,sizeLB,sizeLB});

    Tensor<Scalar,0> E2LRP_ABAB  = E2d_L_evn.contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
    Tensor<Scalar,0> E2LRP_ABBA  = E2d_L_evn.contract(transfer_matrix_LBGA, idx({0,1},{0,1})).contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
    Tensor<Scalar,0> E2LRP_BABA  = E2d_L_odd.contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
    Tensor<Scalar,0> E2LRP_BAAB  = E2d_L_odd.contract(transfer_matrix_LAGB, idx({0,1},{0,1})).contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));


    Scalar e2ab           = E2AB(0);
    Scalar e2ba           = E2BA(0);
    Scalar e2aba_1        = E2ABA_1(0);
    Scalar e2bab_1        = E2BAB_1(0);
    Scalar e2aba_2        = E2ABA_2(0);
    Scalar e2bab_2        = E2BAB_2(0);
    Scalar e2lrpabab      = E2LRP_ABAB(0);
    Scalar e2lrpabba      = E2LRP_ABBA(0);
    Scalar e2lrpbaba      = E2LRP_BABA(0);
    Scalar e2lrpbaab      = E2LRP_BAAB(0);
//    std::cout << setprecision(16);
//    double e_exact = -1.2732395447351625;
//    std::cout << " e_exact          =  "<<  e_exact  << std::endl;
//    std::cout << " e_new            = " <<  e_new           << std::endl;
//    std::cout << " e2_exact         =  "<<  pow(e_exact,2.0) << std::endl;
//    std::cout << " e2ab             = " <<  e2ab            << std::endl;
//    std::cout << " e2ba             = " <<  e2ba            << std::endl;
//    std::cout << " e2aba_1          = " <<  e2aba_1         << std::endl;
//    std::cout << " e2bab_1          = " <<  e2bab_1         << std::endl;
//    std::cout << " e2aba_2          = " <<  e2aba_2         << std::endl;
//    std::cout << " e2bab_2          = " <<  e2bab_2         << std::endl;
//    std::cout << " e2abab_1         = " <<  e2abab_1        << std::endl;
//    std::cout << " e2baba_2         = " <<  e2baba_2        << std::endl;
//    std::cout << " e2baba_3         = " <<  e2baba_3        << std::endl;
//    std::cout << " e2abab_4         = " <<  e2abab_4        << std::endl;
//    std::cout << " e2lrpabab        = " <<  e2lrpabab        << std::endl;
//    std::cout << " e2lrpabba        = " <<  e2lrpabba        << std::endl;
//    std::cout << " e2lrpbaba        = " <<  e2lrpbaba        << std::endl;
//    std::cout << " e2lrpbaab        = " <<  e2lrpbaab        << std::endl;

    Scalar VarE  = 0.5*(e2ab + e2ba) + 0.5*(e2aba_1  + e2bab_1  + e2aba_2  + e2bab_2 )  + e2lrpabab + e2lrpabba + e2lrpbaba  + e2lrpbaab ;
    return real(VarE) ;
}



