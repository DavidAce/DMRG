//
// Created by david on 2017-10-05.
//
//#define EIGEN_USE_MKL_ALL

#include "class_superblock.h"
#include "general/n_math.h"
#include <iomanip>
#include <general/class_svd_wrapper.h>
#include <general/class_eig_spectra_wrapper.h>
using namespace std;
using namespace Textra;

// According to Phien, H. N., Mcculloch, I. P., & Vidal, G. (2015). Fast convergence of imaginary time evolution tensor network algorithms by recycling the environment. Physical Review B - Condensed Matter and Materials Physics, 91(11), 1–18. https://doi.org/10.1103/PhysRevB.91.115137
void class_superblock::canonicalize_iMPS_iterative() {
     if (MPS.L_tail.size() <= 1) {
         return;
     }
     if (MPS.LB.size() != MPS.LA.size()) {
         //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
         return;
     }
     if (MPS.LB.size() < 2 || MPS.LA.size() < 2) {
         //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
         return;
     }
     double diffA = 1;
     double diffB = 1;
     cout << setprecision(16);
     class_SVD<Scalar> svd;
     class_eig<Scalar> eig;
     int iter = 0;
     Tensor<Scalar,3> GA = MPS.GA;
     Tensor<Scalar,3> GB = MPS.GB;
     Tensor<Scalar,1> LA = MPS.LA;
     Tensor<Scalar,1> LB = MPS.LB;
     update_bond_dimensions();

     while (diffA + diffB > 1e-8) {
         iter++;
         //Coarse grain
         Tensor<Scalar,3> GA_old = GA;
         Tensor<Scalar,3> GB_old = GB;
         Tensor<Scalar,1> LA_old = LA;
         Tensor<Scalar,1> LB_old = LB;

         long chiA = LA.size();
         long chiB = LB.size();
         Tensor<Scalar,2> VA_R = GA.contract(asDiagonalSquared(LA) , idx<1>({2},{0}))
                           .contract(GA.conjugate()         , idx<2>({0,2},{0,2}));
         Tensor<Scalar,2> VA_L = GA.contract(asDiagonalSquared(LB) , idx<1>({1},{1}))
                           .contract(GA.conjugate()         , idx<2>({0,2},{0,1}));
         Tensor<Scalar,2> VB_R = GB.contract(asDiagonalSquared(LB) , idx<1>({2},{0}))
                           .contract(GB.conjugate()         , idx<2>({0,2},{0,2}));
         Tensor<Scalar,2> VB_L = GB.contract(asDiagonalSquared(LA) , idx<1>({1},{1}))
                           .contract(GB.conjugate()         , idx<2>({0,2},{0,1}));
         cout << "VA_R: \n" << VA_R << '\n';
         cout << "VA_L: \n" << VA_L << '\n';
         cout << "VB_R: \n" << VB_R << '\n';
         cout << "VB_L: \n" << VB_L << '\n';
         //============ Find square root using Eigenvalue decomposition ============ //
         auto[XA, XAi] = eig.squareRoot_w_inverse(VA_R);
         auto[YA, YAi] = eig.squareRoot_w_inverse(VA_L);
         auto[XB, XBi] = eig.squareRoot_w_inverse(VB_R);
         auto[YB, YBi] = eig.squareRoot_w_inverse(VB_L);

         Tensor<Scalar,2> YA_LA_XB = YA
                 .contract(asDiagonal(LA),idx<1>({1},{0}))
                 .contract(XB,            idx<1>({1},{0}));
         Tensor<Scalar,2> YB_LB_XA = YB
                 .contract(asDiagonal(LB),idx<1>({1},{0}))
                 .contract(XA,            idx<1>({1},{0}));

         Tensor<Scalar,2> UA,VA, UB, VB;
         std::tie(UA,LA,VA) = svd.decompose(YA_LA_XB, chi_max);
         std::tie(UB,LB,VB) = svd.decompose(YB_LB_XA, chi_max);
         GA = VB
                 .contract(XAi,   idx<1>({1},{0}))
                 .contract(GA_old,idx<1>({1},{1}))
                 .contract(YAi,   idx<1>({2},{0}))
                 .contract(UA,    idx<1>({2},{0}))
                 .shuffle(array3{1, 0, 2});
         GB = VA
                 .contract(XBi,    idx<1>({1},{0}))
                 .contract(GB_old, idx<1>({1},{1}))
                 .contract(YBi,    idx<1>({2},{0}))
                 .contract(UB,     idx<1>({2},{0}))
                 .shuffle(array3{1, 0, 2});

         Tensor<Scalar,0> tempA = (LA - LA_old.slice(array1{0}, array1{chiA})).square().sum().sqrt().abs();
         Tensor<Scalar,0> tempB = (LB - LB_old.slice(array1{0}, array1{chiB})).square().sum().sqrt().abs();
         diffA = tempA(0);
         diffB = tempB(0);
     }
     MPS.GA = GA;
     MPS.GB = GB;
     MPS.LA = LA;
     MPS.LB = LB;
     MPS.L_tail = LB;
     cout << "iter: " << iter << endl;

}


//============================================================================//
// Canonicalization, or orthonormalization of the MPS
//============================================================================//


// According to Orus and Vidal
void class_superblock::canonicalize_iMPS(){
    if (MPS.L_tail.size() <= 1) {
        return;
    }
    if (MPS.LB.size() != MPS.LA.size()) {
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    if (MPS.LB.size() < 2 || MPS.LA.size() < 2) {
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    class_SVD<Scalar> svd;
    class_eig<Scalar> eig;

    cout << setprecision(16);
    update_bond_dimensions(); //Gives us chiL __|__chi__|__chiR
    Tensor<Scalar,3> GA = MPS.GA;
    Tensor<Scalar,3> GB = MPS.GB;
    Tensor<Scalar,1> LA = MPS.LA;
    Tensor<Scalar,1> LB = MPS.LB;
    long chiLA = LA.size();
    long chiLB = LB.size();


    Tensor<Scalar,3> GALA = GA.contract(asDiagonal(LA), idx<1>({2},{0}));                          //(eq 107)
    Tensor<Scalar,3> LBGA = asDiagonal(LB).contract(GA, idx<1>({1},{1})).shuffle(array3{1, 0, 2}); //(eq 106)
    Tensor<Scalar,3> GBLB = GB.contract(asDiagonal(LB), idx<1>({2},{0}));                          //(eq 109)
    Tensor<Scalar,3> LAGB = asDiagonal(LA).contract(GB, idx<1>({1},{1})).shuffle(array3{1, 0, 2}); //(eq 108)
    Tensor<Scalar,4> MR = GALA.contract(GALA.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    Tensor<Scalar,4> ML = LBGA.contract(LBGA.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    Tensor<Scalar,4> NR = GBLB.contract(GBLB.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    Tensor<Scalar,4> NL = LAGB.contract(LAGB.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    Tensor<Scalar,2> NRML = NR.contract(ML,idx<2>({2,3},{0,1})).reshape(array2{chiLA*chiLA,chiLA*chiLA});
    Tensor<Scalar,2> NLMR = NL.contract(MR,idx<2>({2,3},{0,1})).reshape(array2{chiLA*chiLA,chiLA*chiLA});
    Tensor<Scalar,2> MRNL = MR.contract(NL,idx<2>({2,3},{0,1})).reshape(array2{chiLB*chiLB,chiLB*chiLB});
    Tensor<Scalar,2> MLNR = ML.contract(NR,idx<2>({2,3},{0,1})).reshape(array2{chiLB*chiLB,chiLB*chiLB});
    cout << "NRML: " << NRML.dimensions() << '\n' << NRML << '\n';
    cout << "NLMR: " << NLMR.dimensions() << '\n' << NLMR << '\n';
    cout << "MRNL: " << MRNL.dimensions() << '\n' << MRNL << '\n';
    cout << "MLNR: " << MLNR.dimensions() << '\n' << MLNR << '\n';
    auto[VA_R, VA_R_val] = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::R>(NRML, 1,
                                                                                                  {chiLA, chiLA});
    auto[VA_L, VA_L_val] = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::L>(NLMR, 1,
                                                                                                  {chiLA, chiLA});
    auto[VB_R, VB_R_val] = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::R>(MRNL, 1,
                                                                                                  {chiLB, chiLB});
    auto[VB_L, VB_L_val] = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::L>(MLNR, 1,
                                                                                                  {chiLB, chiLB});

    cout << "VA_R eigenvalue: " << VA_R_val << "\nVA_R eigenvector: \n" << VA_R << "\n";
    cout << "VA_L eigenvalue: " << VA_L_val << "\nVA_L eigenvector: \n" << VA_L << "\n";
    cout << "VB_R eigenvalue: " << VB_R_val << "\nVB_R eigenvector: \n" << VB_R << "\n";
    cout << "VB_L eigenvalue: " << VB_L_val << "\nVB_L eigenvector: \n" << VB_L << "\n";


    //============ Find square root using Eigenvalue decomposition ============ //
    auto[XA, XAi] = eig.squareRoot_w_inverse(VA_R.real());
    auto[YA, YAi] = eig.squareRoot_w_inverse(VA_L.real());
    auto[XB, XBi] = eig.squareRoot_w_inverse(VB_R.real());
    auto[YB, YBi] = eig.squareRoot_w_inverse(VB_L.real());
    cout << "XA: \n" << XA << '\n';
    cout << "YA: \n" << YA << '\n';
    cout << "XB: \n" << XB << '\n';
    cout << "YB: \n" << YB << '\n';

    Tensor<Scalar,2> YA_LA_XB = YA
            .contract(asDiagonal(LA),idx<1>({1},{0}))
            .contract(XB,            idx<1>({1},{0}));
    Tensor<Scalar,2> YB_LB_XA = YB
            .contract(asDiagonal(LB),idx<1>({1},{0}))
            .contract(XA,            idx<1>({1},{0}));
    Tensor<Scalar,2> UA,VA, UB, VB;
    std::tie(UA,LA,VA) = svd.decompose(YA_LA_XB, chi_max);
    std::tie(UB,LB,VB) = svd.decompose(YB_LB_XA, chi_max);

    cout << "UA: \n" << UA << '\n';
    cout << "LA: \n" << LA << '\n';
    cout << "VA: \n" << VA << '\n';

    GA = VB
            .contract(XAi,   idx<1>({1},{0}))
            .contract(MPS.GA,idx<1>({1},{1}))
            .contract(YAi,   idx<1>({2},{0}))
            .contract(UA,    idx<1>({2},{0}))
            .shuffle(array3{1, 0, 2});
    GB = VA
            .contract(XBi,    idx<1>({1},{0}))
            .contract(MPS.GB, idx<1>({1},{1}))
            .contract(YBi,    idx<1>({2},{0}))
            .contract(UB,     idx<1>({2},{0}))
            .shuffle(array3{1, 0, 2});

    cout << "Check:\n" ;
    GALA = GA.contract(asDiagonal(LA), idx<1>({2},{0}));                          //(eq 107)
    LBGA = asDiagonal(LB).contract(GA, idx<1>({1},{1})).shuffle(array3{1, 0, 2}); //(eq 106)
    GBLB = GB.contract(asDiagonal(LB), idx<1>({2},{0}));                          //(eq 109)
    LAGB = asDiagonal(LA).contract(GB, idx<1>({1},{1})).shuffle(array3{1, 0, 2}); //(eq 108)
    MR   = GALA.contract(GALA.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    ML   = LBGA.contract(LBGA.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    NR   = GBLB.contract(GBLB.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    NL   = LAGB.contract(LAGB.conjugate(),idx<1>({0},{0})).shuffle(array4{0,2,1,3});
    NRML = NR.contract(ML,idx<2>({2,3},{0,1})).reshape(array2{chiLA*chiLA,chiLA*chiLA});
    NLMR = NL.contract(MR,idx<2>({2,3},{0,1})).reshape(array2{chiLA*chiLA,chiLA*chiLA});
    MRNL = MR.contract(NL,idx<2>({2,3},{0,1})).reshape(array2{chiLB*chiLB,chiLB*chiLB});
    MLNR = ML.contract(NR,idx<2>({2,3},{0,1})).reshape(array2{chiLB*chiLB,chiLB*chiLB});
    std::tie(VA_R, VA_R_val) = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::R>(NRML, 1,
                                                                                                      {chiLA, chiLA});
    std::tie(VA_L, VA_L_val) = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::L>(NLMR, 1,
                                                                                                      {chiLA, chiLA});
    std::tie(VB_R, VB_R_val) = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::R>(MRNL, 1,
                                                                                                      {chiLB, chiLB});
    std::tie(VB_L, VB_L_val) = eig.solve_dominant<Spectra::Form::GENERAL, Spectra::Ritz::LARGEST_MAGN, Spectra::Side::L>(MLNR, 1,
                                                                                                      {chiLB, chiLB});

    cout << "VA_R eigenvalue: " << VA_R_val << "\nVA_R eigenvector: \n" << VA_R << "\n";
    cout << "VA_L eigenvalue: " << VA_L_val << "\nVA_L eigenvector: \n" << VA_L << "\n";
    cout << "VB_R eigenvalue: " << VB_R_val << "\nVB_R eigenvector: \n" << VB_R << "\n";
    cout << "VB_L eigenvalue: " << VB_L_val << "\nVB_L eigenvector: \n" << VB_L << "\n";

//    MPS.GA =  GA;
//    MPS.GB =  GB;
//    MPS.LA =  LA;
//    MPS.LB =  LB;
//    MPS.L_tail = LB;
    exit(1);
    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib

}



// According to Orus and Vidal
//void class_superblock::canonicalize2(){
//    if (MPS.L_tail.size() <= 1){
//        return;
//    }
//    class_SVD SVD;
//    SelfAdjointEigenSolver<MatrixXd> esR;
//    SelfAdjointEigenSolver<MatrixXd> esL;
//    //Coarse grain
//    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib
//
//    // G is a Tensor4d  chia__|__|__chib where we expect chiL = chiR for infinite algorithms,
//    // but for finite ones chiL=/=chiR because GA.dimension(1) =/= GB.dimension(2).
//    Tensor4d G = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)})
//            .contract(MPS.GB, idxlist1{idx2(2, 1)}).shuffle(array4{0, 2, 1, 3});
//    Tensor4d G_old = G;
//    //L is simply of size chiR because L = LB;
//    Tensor1d L = MPS.LB;
//
//    Tensor4d GL = G.contract(asDiagonal(L), idxlist1{idx2(3, 0)});
//    Tensor4d LG = asDiagonal(L).contract(G, idxlist1{idx2(1, 2)}).shuffle(array4{1,2,0,3});
//
//    //Watch out for a spurious bug without .eval() here!
//    Tensor2d TR = GL
//            .contract(GL.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiR * chiR, chiR * chiR});
//    Tensor2d TL = LG
//            .contract(LG.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiR * chiR, chiR * chiR});
//    cout << setprecision(16);
//
//
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
//    cout << "TR: \n" << GL.contract(GL.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(3,3)}) << endl;
//    cout << "TL: \n"<<  LG.contract(LG.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(2,2)}) << endl;
//
////    Tensor2d VR = Math::dominant_eigenvector<Math::Handedness::RGHT>(TR).real().reshape(array2{chiR, chiR});
////    Tensor2d VL = Math::dominant_eigenvector<Math::Handedness::RGHT>(TL.shuffle(array2{1,0})).real().reshape(array2{chiR, chiR});
//
//    cout << "VR    =  \n" << VR << '\n';
//    cout << "VL    =  \n" << VL << '\n';
//
//
//    //============================================================================//
//    // Find square root using SVD decomposition
//    //============================================================================//
////    SVD.compute(Tensor2d_to_MatrixXd(VR), Eigen::ComputeThinU | Eigen::ComputeThinV);
////    Tensor2d X               = MatrixXd_to_Tensor2d( SVD.matrixU()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixU().transpose());
////    Tensor2d Xi              = Tensor2d_inverse(X);
////    cout << "X : \n" << X << endl;
////    cout << "Xi : \n" << Xi << endl;
//////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
//////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
////
////    SVD.compute(Tensor2d_to_MatrixXd(VL), Eigen::ComputeThinU | Eigen::ComputeThinV);
////    Tensor2d Y               = MatrixXd_to_Tensor2d( SVD.matrixV().transpose()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixV());
////    Tensor2d Yi              = Tensor2d_inverse(Y);
////    cout << "Y : \n" << Y << endl;
////    cout << "Yi : \n" << Yi << endl;
////        cout << "YTiYT: \n" << YTi.contract(YT, idxlist1{idx2(1,0)}) << endl;
////        cout << "YTY: \n" << YT.contract(Y, idxlist1{idx2(1,0)}) << endl;
//
//
////    ============================================================================//
////     Find square root using Eigenvalue decomposition
////    ============================================================================//
//    esR.compute(Tensor2d_to_MatrixXd(VR));
//    Tensor2d _X    = MatrixXd_to_Tensor2d(esR.operatorSqrt());
//    Tensor2d _XT   = _X.shuffle(array2{1,0});
//    Tensor2d _Xi   = MatrixXd_to_Tensor2d(esR.operatorInverseSqrt());
//    cout << "XXT: \n" << _X.contract(_XT, idxlist1{idx2(1,0)}) << endl;
//
//    esL.compute(Tensor2d_to_MatrixXd(VL));
//    Tensor2d _YT    = MatrixXd_to_Tensor2d(esL.operatorSqrt().transpose());
//    Tensor2d _Y     = MatrixXd_to_Tensor2d(esL.operatorSqrt());
//    Tensor2d _YTi   = MatrixXd_to_Tensor2d(esL.operatorInverseSqrt().transpose());
//    cout << "YTY: \n"   << _YT.contract(_Y, idxlist1{idx2(1,0)}) << endl;
//
//
//    Tensor2d YLX  = _YT.contract(asDiagonal(L), idxlist1{idx2(1,0)})
//            .contract(_X, idxlist1{idx2(1,0)});
//
//
//
////    SVD.compute(Tensor2d_to_MatrixXd(YLX), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    Tensor2d Suzuki_Trotter_2nd_order,V;
//    std::tie(Suzuki_Trotter_2nd_order,L,V) = SVD.decompose(YLX);
////    Tensor2d U = MatrixXd_to_Tensor2d(SVD.matrixU());
////    Tensor2d V = MatrixXd_to_Tensor2d(SVD.matrixV().transpose());
////    L          = MatrixXd_to_Tensor1d(SVD.singularValues().normalized());
//
//    G         = V
//            .contract(_Xi, idxlist1{idx2(1,0)})
//            .contract(G_old, idxlist1{idx2(1,2)})
//            .contract(_YTi, idxlist1{idx2(3, 0)})
//            .contract(Suzuki_Trotter_2nd_order, idxlist1{idx2(3,0)})
//            .shuffle(array4{1,2,0,3});
//
//    cout << "L new: \n" << L << endl;
//
//
//    cout << "\n\n Sanity check \n" << endl;
//    GL = G.contract(asDiagonal(L), idxlist1{idx2(3, 0)});
//    LG = G.contract(asDiagonal(L), idxlist1{idx2(2, 1)}).shuffle(array4{0, 1, 3, 2});
//    //Watch out for a spurious bug without .eval() here!
//    TR = GL
//            .contract(GL.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiR * chiR, chiR * chiR});
//    TL = LG
//            .contract(LG.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiR * chiR, chiR * chiR});
//    cout << "TR    =  \n" << TR << '\n';
//    cout << "TL    =  \n" << TL << '\n';
//    cout << "Transf R: \n" << GL.contract(GL.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(3,3)}) << endl;
//    cout << "Transf L: \n"<<  LG.contract(LG.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(2,2)}) << endl;
//    VR = Math::dominant_eigenvector<Math::Handedness::RGHT>(TR).real().reshape(array2{chiR, chiR});
//    VL = Math::dominant_eigenvector<Math::Handedness::RGHT>(TL.shuffle(array2{1,0})).real().reshape(array2{chiR, chiR});
//
//    cout << "VR    =  \n" << VR << '\n';
//    cout << "VL    =  \n" << VL << '\n';
//
//
//
//
//
//
//    //Separate
//    Tensor2d LGL =   asDiagonal(L)
//            .contract(G, idxlist1{idx2(1,2)})
//            .contract(asDiagonal(L), idxlist1{idx2(3,0)})
//            .shuffle(array4{1,0,2,3})
//            .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});
//
////    SVD.compute(Tensor2d_to_MatrixXd(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    auto[U3, LA, V3] = SVD.schmidt(LGL, d, chiL, chi_max, chiR);
//    truncation_error             = SVD.get_truncation_error();
//
////    chiL                         = std::min(SVD.rank(),chi_max);
////    Tensor3d U3                   = Matrix_to_Tensor<3,double>(SVD.matrixU(),{d,chiR,chiL});
////    MPS.LA                        = Matrix_to_Tensor<1,double>(SVD.singularValues().head(chiL).normalized(), {chiL});
////    Tensor3d V3                   = Matrix_to_Tensor<3,double>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});
//    MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
//    MPS.LA = LA;
//    MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
//    MPS.LB      = Tensor1d_normalize(L);
//    MPS.L_tail  = MPS.LB;
//    cout << "LA new: \n" << MPS.LA << endl;
//    cout << "LB new: \n" << MPS.LB << endl;
////        MPS.swap_AB();
//
//
//    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib
//
//}



// According to Ran, Shi-Ju, Emanuele Tirrito, Cheng Peng, Xi Chen, Gang Su, and Maciej Lewenstein. 2017. “Review of Tensor Network Contraction Approaches,” September. http://arxiv.org/abs/1708.09213.
//void class_superblock::canonicalize4(){
//    if (MPS.L_tail.size() <= 1){
//        return;
//    }
//    cout << setprecision(16);
//
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
//    SelfAdjointEigenSolver<MatrixXd> esR;
//    SelfAdjointEigenSolver<MatrixXd> esL;
//    class_SVD SVD;
//    SVD.setThreshold(settings::precision::SVDThreshold);
//    for(int i = 0; i < 2; i++ ) {
//        update_bond_dimensions(); //Gives us chiL __|__chi__|__chiR
//        //Following the article:
//        //Ran, Shi-Ju, Emanuele Tirrito, Cheng Peng, Xi Chen, Gang Su, and Maciej Lewenstein. 2017. “Review of Tensor Network Contraction Approaches,” September. http://arxiv.org/abs/1708.09213.
//        // Lambda -> Lambda^B
//        // Gamma  -> Lambda^A
//        // A      -> Gamma^A
//        // B      -> Gamma^B
//        long chiA = MPS.LA.size();
//        long chiB = MPS.LB.size();
//        Tensor3d LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//        Tensor3d GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//        Tensor3d LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//        Tensor3d GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//        Tensor4d ML = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//        Tensor4d MR = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//        Tensor4d NL = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//        Tensor4d NR = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//
//        cout << "ML dim: " << ML.dimensions() << endl;
//        cout << "MR dim: " << MR.dimensions() << endl;
//        cout << "NL dim: " << NL.dimensions() << endl;
//        cout << "NR dim: " << NR.dimensions() << endl;
////
////        Tensor2d NRML = NR.contract(ML, idxlist2{idx2(2, 0), idx2(3, 1)})
////                           .shuffle(array4{0,1,2,3})
////                           .reshape(array2{chiA*chiA,chiA*chiA}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)
////        Tensor2d NLMR = NL.contract(MR, idxlist2{idx2(2, 0), idx2(3, 1)})
////                          .shuffle(array4{0,1,2,3})
////                          .reshape(array2{chiA*chiA,chiA*chiA}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)
//
//        Tensor2d NRML = ML.contract(NR, idxlist2{idx2(2, 0), idx2(3, 1)})
//                .shuffle(array4{0,1,2,3})
//                .reshape(array2{chiB*chiB,chiB*chiB}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)
//        Tensor2d NLMR = MR.contract(NL, idxlist2{idx2(2, 0), idx2(3, 1)})
//                .shuffle(array4{0,1,2,3})
//                .reshape(array2{chiB*chiB,chiB*chiB}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)
//
//
////        cout << "NRML: \n" << NRML << '\n';
////        cout << "NLMR: \n" << NLMR << '\n';
//        Tensor2d VR = Math::dominant_eigenvector<Math::Handedness::RGHT>  (NRML).real().reshape(array2{chiB, chiB});
//        Tensor2d VL = Math::dominant_eigenvector<Math::Handedness::RGHT > (NLMR.shuffle(array2{1,0})).real().reshape(array2{chiB, chiB});
//        cout << "VR: \n" << VR << '\n';
//        cout << "VL: \n" << VL << '\n';
//
//        //    ============================================================================//
//        //     Find square root using Eigenvalue decomposition
//        //    ============================================================================//
//        esR.compute(Tensor2d_to_MatrixXd(VR));
//        Tensor2d X    = MatrixXd_to_Tensor2d(esR.operatorSqrt());
//        Tensor2d XT   = X.shuffle(array2{1,0});
//        Tensor2d Xi   = MatrixXd_to_Tensor2d(esR.operatorInverseSqrt());
//        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//        esL.compute(Tensor2d_to_MatrixXd(VL));
//        Tensor2d YT    = MatrixXd_to_Tensor2d(esL.operatorSqrt().transpose());
//        Tensor2d Y     = MatrixXd_to_Tensor2d(esL.operatorSqrt());
//        Tensor2d YTi   = MatrixXd_to_Tensor2d(esL.operatorInverseSqrt().transpose());
//        cout << "YTY: \n"   << YT.contract(Y, idxlist1{idx2(1,0)}) << endl;
//
//
//        cout << "X      =  \n" << X << '\n';
//        cout << "X*XT   =  \n" << X.contract(XT, idxlist1{idx2(1, 0)}) << '\n';
//        cout << "Xi*X   =  \n" << Xi.contract(X, idxlist1{idx2(1, 0)}) << '\n';
//        cout << "Y      =  \n" << Y << '\n';
//        cout << "YT*Y   =  \n" << YT.contract(Y, idxlist1{idx2(1, 0)}) << '\n';
//        cout << "YT*YTi   =  \n" << YTi.contract(YT, idxlist1{idx2(1, 0)}) << '\n';
//
//
//        Tensor2d X_LA_Y = X
//                .contract(asDiagonal(MPS.LA), idxlist1{idx2(1, 0)})
//                .contract(Y, idxlist1{idx2(1, 0)});
//
//        auto[UA,LA,VA] = SVD.decompose(X_LA_Y);
//        MPS.LA = LA;
//
//
//        Tensor3d GA = MPS.GA.contract(Xi, idxlist1{idx2(2, 0)})
//                .contract(UA, idxlist1{idx2(2, 0)});
//        Tensor3d GB = MPS.GB.contract(YTi, idxlist1{idx2(1, 1)})
//                .contract(VA, idxlist1{idx2(2, 1)})
//                .shuffle(array3{0,2,1});
//        MPS.GA = GA;
//        MPS.GB = GB;
//        cout << "LA new: \n" << MPS.LA << endl;
//        cout << "LB new: \n" << MPS.LB << endl;
//        MPS.L_tail = MPS.LB;
//
//
//        cout << "\n\n Sanity Check " << endl;
//        LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//        GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//        LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//        GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//        Tensor2d ML_ = LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1,1)});
//        Tensor2d MR_ = GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2,2)});
//        Tensor2d NL_ = LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1,1)});
//        Tensor2d NR_ = GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2,2)});
//        cout << "ML: \n " <<ML_<< '\n';
//        cout << "MR: \n " <<MR_<< '\n';
//        cout << "NL: \n " <<NL_<< '\n';
//        cout << "MR: \n " <<MR_<< '\n';
//
//
//
//
//
//
//        MPS.swap_AB();
//    }
//
//
//
//    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib
//}


// According to Phien, H. N., Mcculloch, I. P., & Vidal, G. (2015). Fast convergence of imaginary time evolution tensor network algorithms by recycling the environment. Physical Review B - Condensed Matter and Materials Physics, 91(11), 1–18. https://doi.org/10.1103/PhysRevB.91.115137
//void class_superblock::canonicalize3(){
//    if (MPS.L_tail.size() <= 1){
//        return;
//    }
//
//    double diff = 1;
//    cout << setprecision(16) << endl;
//    cout << "LA before: \n" << MPS.LA << endl;
//    cout << "LB before: \n" << MPS.LB << endl;
//    Eigen::BDCSVD<MatrixXd> SVD;
//    SelfAdjointEigenSolver<MatrixXd> esR;
//    SelfAdjointEigenSolver<MatrixXd> esL;
//
//    SVD.setThreshold(1e-16);
//    int iter = 0;
//    while(diff > 1e-8) {
////        cout << "diff: " << diff << endl;
//        iter++;
//        //Coarse grain
//        Tensor4d G  =  MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2,0)})
//                .contract(MPS.GB, idxlist1{idx2(2,1)}).shuffle(array4{0,2,1,3});
//        Tensor1d L       = MPS.LB;
//
//
//        Tensor1d L_old = L;
//        Tensor4d G_old = G;
//
//        Tensor4d GL = G.contract(asDiagonal(L), idxlist1{idx2(3, 0)});
//        Tensor4d LG = asDiagonal(L).contract(G, idxlist1{idx2(1, 2)}).shuffle(array4{1,2,0,3});
//
//        Tensor2d VR = GL.contract(GL.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(3, 3)});
//        Tensor2d VL = LG.contract(LG.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(2, 2)});
//
//        cout << "VR    =  \n" << VR << '\n';
//        cout << "VL    =  \n" << VL << '\n';
//        //============================================================================//
//        // Find square root using SVD decomposition
//        //============================================================================//
//        cout << "\n\nSVD:\n";
//        SVD.compute(Tensor2_to_Matrix<double>(VR), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        Tensor2d X               = matrix_to_tensor2<double>( SVD.matrixU()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixU().transpose());
//        Tensor2d XT              = X.shuffle(array2{1,0});
//        Tensor2d Xi              = matrix_to_tensor2<double>(Tensor2_to_Matrix(X).inverse());
//        cout << "X : \n" << X << endl;
//        cout << "Xi : \n" << Xi << endl;
////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//        SVD.compute(Tensor2_to_Matrix(VL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        Tensor2d Y               = matrix_to_tensor2<double>( SVD.matrixU()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixU().transpose());
//        Tensor2d Yi              = matrix_to_tensor2<double>(Tensor2_to_Matrix(Y).inverse());
//        cout << "Y : \n" << Y << endl;
//        cout << "Yi : \n" << Yi << endl;
////        cout << "YTiYT: \n" << YTi.contract(YT, idxlist1{idx2(1,0)}) << endl;
////        cout << "YTY: \n" << YT.contract(Y, idxlist1{idx2(1,0)}) << endl;
//
//
//        //============================================================================//
//        // Find square root using Eigenvalue decomposition
//        //============================================================================//
//        cout << "\n\nEigenvalue:\n";
//        esR.compute(Tensor2_to_Matrix(VR));
//        Tensor2d _X    = matrix_to_tensor2<double>( esR.eigenvectors()*esR.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt());
//        Tensor2d _XT   = matrix_to_tensor2<double>( esR.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt() * esR.eigenvectors().transpose());
//        Tensor2d _Xi   = matrix_to_tensor2<double>(Tensor2_to_Matrix(X).inverse());
//        cout << "X : \n" << _X << endl;
//        cout << "Xi : \n" << _Xi << endl;
////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//        esL.compute(Tensor2_to_Matrix(VL));
//        Tensor2d _Y              = matrix_to_tensor2<double>( esL.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt() * esL.eigenvectors().transpose());
//        Tensor2d _YT             = matrix_to_tensor2<double>( esL.eigenvectors()*esL.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt());
//        Tensor2d _YTi            = matrix_to_tensor2<double>(Tensor2_to_Matrix(_YT).inverse());
//        cout << "YT : \n"  << _YT << endl;
//        cout << "YTi : \n" << _YTi << endl;
////        cout << "YTiYT: \n" << YTi.contract(YT, idxlist1{idx2(1,0)}) << endl;
////        cout << "YTY: \n"   << YT.contract(Y, idxlist1{idx2(1,0)}) << endl;
//
//
//        Tensor2d YLX = Y.contract(asDiagonal(L), idxlist1{idx2(1,0)})
//                .contract(X, idxlist1{idx2(1,0)});
//
//
//
//        SVD.compute(Tensor2_to_Matrix(YLX), Eigen::ComputeThinU | Eigen::ComputeThinV);
//
//        Tensor2d Suzuki_Trotter_2nd_order = matrix_to_tensor2<double>(SVD.matrixU());
//        Tensor2d V = matrix_to_tensor2<double>(SVD.matrixV().transpose());
//        L         = Matrix_to_Tensor1<double>(SVD.singularValues().normalized());
//
//        G         = V
//                   .contract(Xi, idxlist1{idx2(1,0)})
//                   .contract(G_old, idxlist1{idx2(1,2)})
//                   .contract(Yi, idxlist1{idx2(3, 0)})
//                   .contract(Suzuki_Trotter_2nd_order, idxlist1{idx2(3,0)})
//                   .shuffle(array4{1,2,0,3});
//
//
//        //Separate
//        Tensor2d LGL =   asDiagonal(L)
//                .contract(G, idxlist1{idx2(1,2)})
//                .contract(asDiagonal(L), idxlist1{idx2(3,0)})
//                .shuffle(array4{1,0,2,3})
//                .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});
//
//        SVD.compute(Tensor2_to_Matrix(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        chiL                         = std::min(SVD.rank(),chi_max);
//        truncation_error             = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiL).sum();
//        Tensor3d U3                   = matrix_to_tensor<3,double>(SVD.matrixU(),{d,chiR,chiL});
//        MPS.LA                        = matrix_to_tensor<1,double>(SVD.singularValues().head(chiL).normalized(), {chiL});
//        Tensor3d V3                   = matrix_to_tensor<3,double>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});
//
//        MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
//        MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
//        MPS.LB = L;
//        MPS.L_tail = L;
//
//
//
//        Tensor0d temp = (L-L_old.slice(array1{0}, array1{chiR})).square().sum().sqrt();
//        diff = temp(0);
//        cout << "L: \n" << L << endl;
//        cout << "L diff: " << diff << endl;
//        if(iter > 100 ){
//            exit(0);
//        }
//    }
//
////    Tensor2d LGL =   asDiagonal(L)
////                   .contract(G, idxlist1{idx2(1,2)})
////                   .contract(asDiagonal(L), idxlist1{idx2(3,0)})
////                   .shuffle(array4{1,0,2,3})
////                   .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});
////
////
//////    //Note the difference! U3 is a d x [chib] x chi matrix! Previously we had [chia] in normal SVD truncation
////    SVD.compute(tensor2_to_matrix(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
////    chiL                         = std::min(SVD.rank(),chi_max);
////    truncation_error             = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiL).sum();
////    Tensor3d U3                   = matrix_to_tensor<3>(SVD.matrixU(),{d,chiR,chiL});
////    MPS.LA                       = matrix_to_tensor<1>(SVD.singularValues().head(chiL).normalized(), {chiL});
////    Tensor3d V3                   = matrix_to_tensor<3>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});
//////
////    MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
////    MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
////    MPS.LB = L;
////    MPS.L_tail = L;
//
//    cout << "\n\niter: " << iter << endl;
////    cout << "Full LA after: \n" << SVD.singularValues().normalized() << endl;
//    cout << "LA after: \n" << MPS.LA << endl;
//    cout << "LB after: \n" << MPS.LB << endl;
//
//
//
//    return;
//
//}




//void class_superblock::transfer_matrix_eigenvalues() {
//
//    if (MPS.LB.size()!= MPS.LA.size()){
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    if (MPS.LB.size() < 2 || MPS.LA.size() < 2 ){
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    cout << setprecision(16);
//    update_bond_dimensions();
//    long chiA = MPS.LA.size();
//    long chiB = MPS.LB.size();
//    Tensor3d GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//    Tensor3d GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//    Tensor3d LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//    Tensor3d LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//
//    Tensor4d TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//
////    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
////    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
////    Tensor2d VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
////    Tensor2d VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
//    Spectra::DenseGenMatProd<double> op_TA_R(Tensor_to_Matrix<4>(TA_R, chiB*chiB,chiA*chiA));
//    Spectra::DenseGenMatProd<double> op_TA_L(Tensor_to_Matrix<4>(TA_L, chiB*chiB,chiA*chiA));
//    Spectra::DenseGenMatProd<double> op_TB_R(Tensor_to_Matrix<4>(TB_R, chiA*chiA,chiB*chiB));
//    Spectra::DenseGenMatProd<double> op_TB_L(Tensor_to_Matrix<4>(TB_L, chiA*chiA,chiB*chiB));
//    int ncv = std::min(settings::precision::Spectra::max_ncv,4);
//    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTA_R(&op_TA_R, 2, ncv);
//    eigsTA_R.init();
//    eigsTA_R.compute(10000, 1e-12, Spectra::LARGEST_REAL);
//    std::cout << "Eigenvalue TA_R:\n" << eigsTA_R.eigenvalues()<< std::endl;
//    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTA_L(&op_TA_L, 2, ncv);
//    eigsTA_L.init();
//    eigsTA_L.compute(10000, 1e-12, Spectra::LARGEST_REAL);
//    std::cout << "Eigenvalue TA_L:\n" << eigsTA_L.eigenvalues()<< std::endl;
//    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTB_R(&op_TB_R, 2, ncv);
//    eigsTB_R.init();
//    eigsTB_R.compute(10000, 1e-12, Spectra::LARGEST_REAL);
//    std::cout << "Eigenvalue TB_R:\n" << eigsTB_R.eigenvalues()<< std::endl;
//    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTB_L(&op_TB_L, 2, ncv);
//    eigsTB_L.init();
//    eigsTB_L.compute(10000, 1e-12, Spectra::LARGEST_REAL);
//    std::cout << "Eigenvalue TB_L: \n" << eigsTB_L.eigenvalues()<< std::endl;
//
//    return;
//}

//void class_superblock::canonicalize_iMPS_iterative_again() {
//    if (MPS.L_tail.size() <= 1) {
//        return;
//    }
//    if (MPS.LB.size() != MPS.LA.size()) {
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    if (MPS.LB.size() < 2 || MPS.LA.size() < 2) {
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    double diffA = 1;
//    double diffB = 1;
//    cout << setprecision(16) << endl;
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
//    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//
//    Eigen::BDCSVD<MatrixXd> SVD;
////    SelfAdjointEigenSolver<MatrixXd> es;
////
//    SVD.setThreshold(1e-24);
//    int iter = 0;
//    Tensor3d GA = MPS.GA;
//    Tensor3d GB = MPS.GB;
//    Tensor1d LA = MPS.LA;
//    Tensor1d LB = MPS.LB;
//    update_bond_dimensions();
//
//    while (diffA + diffB > 1e-10) {
////        cout << "diff: " << diff << endl;
//        iter++;
//        //Coarse grain
//        Tensor3d GA_old = GA;
//        Tensor3d GB_old = GB;
//        Tensor1d LA_old = LA;
//        Tensor1d LB_old = LB;
//
//        long chiA = LA.size();
//        long chiB = LB.size();
//        Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//        Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//        Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//        Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//
//        Tensor2d VA_R = GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)});
//        Tensor2d VA_L = LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)});
//        Tensor2d VB_R = GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)});
//        Tensor2d VB_L = LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)});
//
//        //============================================================================//
//        // Find square root using Eigenvalue decomposition
//        //============================================================================//
//        SelfAdjointEigenSolver<MatrixXd> es;
//        es.compute(Tensor2d_to_MatrixXd(VA_R));
//        Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
//        Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//        es.compute(Tensor2d_to_MatrixXd(VA_L));
//        Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//        Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//        es.compute(Tensor2d_to_MatrixXd(VB_R));
//        Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//        Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//        es.compute(Tensor2d_to_MatrixXd(VB_L));
//        Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//        Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//
//        //Symmetric
//        Tensor2d YA_LA_XB = YA
//                .contract(asDiagonal(LA), idxlist1{idx2(1, 0)})
//                .contract(XB, idxlist1{idx2(1, 0)});
//        //Symmetric
//        Tensor2d YB_LB_XA = YB
//                .contract(asDiagonal(LB), idxlist1{idx2(1, 0)})
//                .contract(XA, idxlist1{idx2(1, 0)});
//
//        SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        chiA = std::min(SVD.rank(), chi_max);
//        Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
//        Tensor2d VA = MatrixXd_to_Tensor2d(
//                SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, Suzuki_Trotter_2nd_order and V are not symmetric!
//        LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());
//
//        SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        chiB = std::min(SVD.rank(), chi_max);
//        Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
//        Tensor2d VB = MatrixXd_to_Tensor2d(
//                SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, Suzuki_Trotter_2nd_order and V are not symmetric!
//        LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());
//
////    cout << "UA: \n "<< UA << endl;
////    cout << "VA: \n "<< VA << endl;
////    cout << "UB: \n "<< UB << endl;
////    cout << "VB: \n "<< VB << endl;
//        GA = VB
//                .contract(XAi, idxlist1{idx2(1, 0)})
//                .contract(GA_old, idxlist1{idx2(1, 1)})
//                .contract(YAi, idxlist1{idx2(2, 0)})
//                .contract(UA, idxlist1{idx2(2, 0)})
//                .shuffle(array3{1, 0, 2});
//        GB = VA
//                .contract(XBi, idxlist1{idx2(1, 0)})
//                .contract(GB_old, idxlist1{idx2(1, 1)})
//                .contract(YBi, idxlist1{idx2(2, 0)})
//                .contract(UB, idxlist1{idx2(2, 0)})
//                .shuffle(array3{1, 0, 2});
//
//
//        Tensor0d tempA = (LA - LA_old.slice(array1{0}, array1{chiA})).square().sum().sqrt();
//        Tensor0d tempB = (LB - LB_old.slice(array1{0}, array1{chiB})).square().sum().sqrt();
//        diffA = tempA(0);
//        diffB = tempB(0);
////         cout << "LA diff: " << diffA << endl;
////         cout << "LB diff: " << diffB << endl;
////        if(iter > 100 ){
////            exit(0);
////        }
//    }
//    Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//    Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//    Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//    Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//    cout << "TA_R: \n " << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
//    cout << "TA_L: \n " << LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
//    cout << "TB_R: \n " << GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
//    cout << "TB_L: \n " << LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
////    Tensor2d LGL =   asDiagonal(L)
////                   .contract(G, idxlist1{idx2(1,2)})
////                   .contract(asDiagonal(L), idxlist1{idx2(3,0)})
////                   .shuffle(array4{1,0,2,3})
////                   .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});
////
////
//////    //Note the difference! U3 is a d x [chib] x chi matrix! Previously we had [chia] in normal SVD truncation
////    SVD.compute(tensor2_to_matrix(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
////    chiL                         = std::min(SVD.rank(),chi_max);
////    truncation_error             = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiL).sum();
////    Tensor3d U3                   = matrix_to_tensor<3>(SVD.matrixU(),{d,chiR,chiL});
////    MPS.LA                       = matrix_to_tensor<1>(SVD.singularValues().head(chiL).normalized(), {chiL});
////    Tensor3d V3                   = matrix_to_tensor<3>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});
//////
////    MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
////    MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
////    MPS.LB = L;
////    MPS.L_tail = L;
//
////    cout << "Full LA after: \n" << SVD.singularValues().normalized() << endl;
//    cout << "LA after: \n" << LA << endl;
//    cout << "LB after: \n" << LB << endl;
//    cout << "Entropy A: " << -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
//    cout << "Entropy B: " << -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
//    cout << "\n\niter: " << iter << endl;
//
////     MPS.GA = GA;
////     MPS.GB = GB;
////     MPS.LA = LA;
////     MPS.LB = LB;
////     MPS.L_tail = LB;
////    TRY USING ROBUST SVD!
////             Stoudenmire, E. M., & White, S. R. (2013). Real-space parallel density matrix renormalization group. Physical Review B, 87(April), 155137. https://doi.org/10.1103/PhysRevB.87.155137
//
//}

//void class_superblock::canonicalize_iMPS_vidal() {
//    if (MPS.L_tail.size() <= 1) {
//        return;
//    }
//    if (MPS.LB.size() != MPS.LA.size()) {
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    if (MPS.LB.size() < 2 || MPS.LA.size() < 2) {
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
////    double diffA = 1;
////    double diffB = 1;
//    cout << setprecision(16) << endl;
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
//    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//
//    Eigen::BDCSVD<MatrixXd> SVD;
////    SelfAdjointEigenSolver<MatrixXd> es;
////
//    SVD.setThreshold(1e-24);
//    int iter = 0;
//    Tensor3d GA = MPS.GA;
//    Tensor3d GB = MPS.GB;
//    Tensor1d LA = MPS.LA;
//    Tensor1d LB = MPS.LB;
//    update_bond_dimensions();
//
//    Tensor3d GA_old = GA;
//    Tensor3d GB_old = GB;
//    Tensor1d LA_old = LA;
//    Tensor1d LB_old = LB;
//
//    long chiA = LA.size();
//    long chiB = LB.size();
//    Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//    Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//    Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//    Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//
//    Tensor4d TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//
//
//    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
//    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
//    Tensor2d VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
//    Tensor2d VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
//
//    VA_R = VA_R/VA_R.coeff(0);
//    VA_L = VA_L/VA_L.coeff(0);
//    VB_R = VB_R/VB_R.coeff(0);
//    VB_L = VB_L/VB_L.coeff(0);
//    cout << "VA_R: \n " << VA_R << '\n';
//    cout << "VA_L: \n " << VA_L << '\n';
//    cout << "VB_R: \n " << VB_R << '\n';
//    cout << "VB_L: \n " << VB_L << '\n';
//    //============================================================================//
//    // Find square root using Eigenvalue decomposition
//    //============================================================================//
//    SelfAdjointEigenSolver<MatrixXd> es;
//    es.compute(Tensor2d_to_MatrixXd(VA_R));
//    Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VA_L));
//    Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VB_R));
//    Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VB_L));
//    Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//
//    //Symmetric
//    Tensor2d YA_LA_XB = YA
//            .contract(asDiagonal(LA), idxlist1{idx2(1, 0)})
//            .contract(XB, idxlist1{idx2(1, 0)});
//    //Symmetric
//    Tensor2d YB_LB_XA = YB
//            .contract(asDiagonal(LB), idxlist1{idx2(1, 0)})
//            .contract(XA, idxlist1{idx2(1, 0)});
//
//    SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    chiA = std::min(SVD.rank(), chi_max);
//    Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
//    Tensor2d VA = MatrixXd_to_Tensor2d(
//            SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, Suzuki_Trotter_2nd_order and V are not symmetric!
//    LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());
//
//    SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    chiB = std::min(SVD.rank(), chi_max);
//    Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
//    Tensor2d VB = MatrixXd_to_Tensor2d(
//            SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, Suzuki_Trotter_2nd_order and V are not symmetric!
//    LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());
//
////    cout << "UA: \n "<< UA << endl;
////    cout << "VA: \n "<< VA << endl;
////    cout << "UB: \n "<< UB << endl;
////    cout << "VB: \n "<< VB << endl;
//    GA = VB
//            .contract(XAi, idxlist1{idx2(1, 0)})
//            .contract(GA_old, idxlist1{idx2(1, 1)})
//            .contract(YAi, idxlist1{idx2(2, 0)})
//            .contract(UA, idxlist1{idx2(2, 0)})
//            .shuffle(array3{1, 0, 2});
//    GB = VA
//            .contract(XBi, idxlist1{idx2(1, 0)})
//            .contract(GB_old, idxlist1{idx2(1, 1)})
//            .contract(YBi, idxlist1{idx2(2, 0)})
//            .contract(UB, idxlist1{idx2(2, 0)})
//            .shuffle(array3{1, 0, 2});
//
//
//    Tensor0d tempA = (LA - LA_old.slice(array1{0}, array1{chiA})).square().sum().sqrt();
//    Tensor0d tempB = (LB - LB_old.slice(array1{0}, array1{chiB})).square().sum().sqrt();
//
////         cout << "LA diff: " << diffA << endl;
////         cout << "LB diff: " << diffB << endl;
////        if(iter > 100 ){
////            exit(0);
////        }
//
//
//    cout << "LA after: \n" << LA << endl;
//    cout << "LB after: \n" << LB << endl;
//    cout << "Entropy A: " << -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
//    cout << "Entropy B: " << -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
//
//
//    GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//    LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//    GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//    LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//
//    TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//
//    VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
//    VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
//    VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
//    VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
//    cout << "VA_R: \n " << VA_R << '\n';
//    cout << "VA_L: \n " << VA_L << '\n';
//    cout << "VB_R: \n " << VB_R << '\n';
//    cout << "VB_L: \n " << VB_L << '\n';
//
//    cout << "TA_R: \n " << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
//    cout << "TA_L: \n " << LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
//    cout << "TB_R: \n " << GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
//    cout << "TB_L: \n " << LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
//    cout << "\n\niter: " << iter << endl;
//
//     MPS.GA = GA;
//     MPS.GB = GB;
//     MPS.LA = LA;
//     MPS.LB = LB;
//     MPS.L_tail = LB;
////    TRY USING ROBUST SVD!
////             Stoudenmire, E. M., & White, S. R. (2013). Real-space parallel density matrix renormalization group. Physical Review B, 87(April), 155137. https://doi.org/10.1103/PhysRevB.87.155137
//
//}


//void class_superblock::canonicalize_iMPS_mcculloch() {
//    //Start by creating a 2-site overlap
//
//    if (MPS.LB.size()!= MPS.LA.size()){
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    if (MPS.LB.size() < 2 || MPS.LA.size() < 2 ){
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    Eigen::BDCSVD<MatrixXd> SVD;
//    SVD.setThreshold(1e-14);
//    cout << setprecision(16);
//    update_bond_dimensions();
//    long chiA = MPS.LA.size();
//    long chiB = MPS.LB.size();
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
//    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//
//    Tensor2d P      = asDiagonal(MPS.LB).contract(asInverseDiagonal(MPS.LB), idxlist1{idx2(1,0)});
//    Tensor2d Q      = asInverseDiagonal(MPS.LB).contract(asDiagonal(MPS.LB), idxlist1{idx2(1,0)});
//    Tensor4d thetaP =  asDiagonal(MPS.LB)
//                      .contract(MPS.GA, idxlist1{idx2(1,1)})
//                      .contract(asDiagonal(MPS.LA), idxlist1{idx2(2,0)})
//                      .contract(MPS.GB, idxlist1{idx2(2,1)})
//                      .contract(P, idxlist1{idx2(3,0)});
//    Tensor4d Qtheta =  Q
//            .contract(MPS.GA, idxlist1{idx2(1,1)})
//            .contract(asDiagonal(MPS.LA), idxlist1{idx2(2,0)})
//            .contract(MPS.GB, idxlist1{idx2(2,1)})
//            .contract(asDiagonal(MPS.LB), idxlist1{idx2(3,0)});
//
//    Tensor4d TL     = thetaP.contract(thetaP,idxlist2{idx2(1,1),idx2(2,2)}).shuffle(array4{0,2,1,3});
//    Tensor4d TR     = Qtheta.contract(Qtheta,idxlist2{idx2(1,1),idx2(2,2)}).shuffle(array4{0,2,1,3});
//
//    Tensor2d VL     = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TL.reshape(array2{chiB*chiB,chiB*chiB})).real().reshape(array2{chiB, chiB});
//    Tensor2d VR     = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TR.reshape(array2{chiA*chiA,chiA*chiA})).real().reshape(array2{chiA, chiA});
//
//    //Rescale
//
//    cout << "VL : \n" << VL << endl;
//    cout << "VR : \n" << VR << endl;
//    SelfAdjointEigenSolver<MatrixXd> es;
//    es.compute(Tensor2d_to_MatrixXd(VL));
//    Tensor2d X  = MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d Xi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VR));
//    Tensor2d Y  = MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d Yi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    cout << "X : \n" << X << endl;
//    cout << "Y : \n" << Y << endl;
//
//    Tensor2d EL = Xi.contract(VL,idxlist1{idx2(1,0)}).contract(Xi,idxlist1{idx2(1,0)});
//    cout << "EL: \n" << EL << endl;
//
//    Tensor2d LB = X.contract(asDiagonal(MPS.LB), idxlist1{idx2(1,0)}).contract(Y,idxlist1{idx2(1,0)});
//    cout << "LB: \n" << LB << endl;
//    cout << "S_E: " <<  -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//
//
//}


//void class_superblock::canonicalize_iMPS() {
//    if (MPS.LB.size()!= MPS.LA.size()){
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    if (MPS.LB.size() < 2 || MPS.LA.size() < 2 ){
//        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
//        return;
//    }
//    cout << setprecision(16);
//    update_bond_dimensions();
//    long chiA = MPS.LA.size();
//    long chiB = MPS.LB.size();
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
//    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//
//    Eigen::BDCSVD<MatrixXd> SVD;
//    SVD.setThreshold(1e-14);
//    Tensor3d GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//    Tensor3d GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//    Tensor3d LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//    Tensor3d LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//
//    Tensor4d TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    Tensor4d TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
////
////    Tensor4d MRNR = TA_R.contract(TB_R.conjugate(), idxlist2{idx2(2, 0), idx2(3,1)});
////    Tensor4d MLNL = TA_L.contract(TB_L.conjugate(), idxlist2{idx2(2, 2), idx2(1,3)});
//////    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RIGHT>  (MRNR.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
//////    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::RIGHT> (MLNL.reshape(array2{chiB*chiB,chiA*chiA}).shuffle(array2{1,0})).real().reshape(array2{chiA, chiA});
//////
//
//    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
//    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
//    Tensor2d VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
//    Tensor2d VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
//
//    cout << "VA_R : \n" << VA_R << endl;
//    cout << "VA_L : \n" << VA_L << endl;
//    cout << "VB_R : \n" << VB_R << endl;
//    cout << "VB_L : \n" << VB_L << endl;
//
////    SVD.compute(Tensor2d_to_MatrixXd(VA_R), Eigen::ComputeThinU);
////    Tensor2d XA  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    Tensor2d XAi = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    SVD.compute(Tensor2d_to_MatrixXd(VA_L), Eigen::ComputeThinU);
////    Tensor2d YA   = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    Tensor2d YAi  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    SVD.compute(Tensor2d_to_MatrixXd(VB_R), Eigen::ComputeThinU);
////    Tensor2d XB   = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    Tensor2d XBi  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    SVD.compute(Tensor2d_to_MatrixXd(VB_L), Eigen::ComputeThinU);
////    Tensor2d YB   = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////    Tensor2d YBi  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
////
//
//    //============================================================================//
//    // Find square root using Eigenvalue decomposition
//    //============================================================================//
//    SelfAdjointEigenSolver<MatrixXd> es;
//    es.compute(Tensor2d_to_MatrixXd(VA_R));
//    Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VA_L));
//    Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VB_R));
//    Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//    es.compute(Tensor2d_to_MatrixXd(VB_L));
//    Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
//    Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
//
////    cout << "XA         =  \n"  << XA << '\n';
////    cout << "XA*XAi     =  \n"  << XA.contract(XAi, idxlist1{idx2(1, 0)}) << '\n';
////    cout << "YA         =  \n"  << YA << '\n';
////    cout << "YAi*YA     =  \n"  << YAi.contract(YA, idxlist1{idx2(1, 0)}) << '\n';
////
////    cout << "XB         =  \n"  << XB << '\n';
////    cout << "XB*XBi     =  \n"  << XB.contract(XBi, idxlist1{idx2(1, 0)}) << '\n';
////    cout << "YB         =  \n"  << YB << '\n';
////    cout << "YBi*YB     =  \n"  << YBi.contract(YB, idxlist1{idx2(1, 0)}) << '\n';
////
////    cout << "XA dims: " << XA.dimensions() << endl;
////    cout << "YA dims: " << YA.dimensions() << endl;
////    cout << "XB dims: " << XB.dimensions() << endl;
////    cout << "YB dims: " << YB.dimensions() << endl;
//
//    //Symmetric
//    Tensor2d YA_LA_XB = YA
//            .contract(asDiagonal(MPS.LA), idxlist1{idx2(1, 0)})
//            .contract(XB, idxlist1{idx2(1, 0)});
//    //Symmetric
//    Tensor2d YB_LB_XA = YB
//            .contract(asDiagonal(MPS.LB), idxlist1{idx2(1, 0)})
//            .contract(XA, idxlist1{idx2(1, 0)});
//
//    SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    chiA        = std::min(SVD.rank(),chi_max);
//    Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
//    Tensor2d VA = MatrixXd_to_Tensor2d(SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, Suzuki_Trotter_2nd_order and V are not symmetric!
//    Tensor1d LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());
//
//    SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    chiB        = std::min(SVD.rank(),chi_max);
//    Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
//    Tensor2d VB = MatrixXd_to_Tensor2d(SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, Suzuki_Trotter_2nd_order and V are not symmetric!
//    Tensor1d LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());
//    cout << "LA new: \n" << LA << endl;
//    cout << "LB new: \n" << LB << endl;
//    cout << "S_E A: " <<  -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//    cout << "S_E B: " <<  -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
//
////    cout << "UA: \n "<< UA << endl;
////    cout << "VA: \n "<< VA << endl;
////    cout << "UB: \n "<< UB << endl;
////    cout << "VB: \n "<< VB << endl;
//    Tensor3d GA = VB
//               .contract(XAi, idxlist1{idx2(1,0)})
//               .contract(MPS.GA, idxlist1{idx2(1,1)})
//               .contract(YAi, idxlist1{idx2(2, 0)})
//               .contract(UA, idxlist1{idx2(2,0)})
//               .shuffle(array3{1,0,2});
//    Tensor3d GB = VA
//            .contract(XBi, idxlist1{idx2(1,0)})
//            .contract(MPS.GB, idxlist1{idx2(1,1)})
//            .contract(YBi, idxlist1{idx2(2, 0)})
//            .contract(UB, idxlist1{idx2(2,0)})
//            .shuffle(array3{1,0,2});
//
//    GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
//    LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
//    GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
//    LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
//    TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//    TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//
//
//    VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
//    VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
//    VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
//    VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
//
//
//    cout << "VA_R: \n " << VA_R << '\n';
//    cout << "VA_L: \n " << VA_L << '\n';
//    cout << "VB_R: \n " << VB_R << '\n';
//    cout << "VB_L: \n " << VB_L << '\n';
//    //Check full 2-site Transfer matrix
//    MatrixXd IL = MatrixXd::Identity(TA_L.dimension(0), TA_L.dimension(1));
//    MatrixXd IR = MatrixXd::Identity(TB_R.dimension(2), TB_R.dimension(3));
//    Tensor2d TL = TA_L.contract(TB_L, idxlist2{idx2(2,0), idx2(3,1)}).contract(MatrixXd_to_Tensor2d(IL), idxlist2{idx2(0,0),idx2(1,1)});
//    Tensor2d TR = TA_R.contract(TB_R, idxlist2{idx2(2,0), idx2(3,1)}).contract(MatrixXd_to_Tensor2d(IR), idxlist2{idx2(2,0),idx2(3,1)});
//
//    cout << "TL: \n" << TL << endl;
//    cout << "TR: \n" << TR << endl;
////    MPS.GA = GA;
////    MPS.GB = GB;
////    MPS.LA = LA;
////    MPS.LB = LB;
////    MPS.L_tail = LB;
//
//}


