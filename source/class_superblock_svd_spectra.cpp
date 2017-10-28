//
// Created by david on 8/15/17.
//

#include <class_superblock.h>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <SymEigsSolver.h>
#include <GenEigsSolver.h>
#include <iomanip>
#include <n_settings.h>
#include <n_math.h>

using namespace std;
using namespace Textra;
using namespace Eigen;


//============================================================================//
// Find smallest eigenvalue using Spectra.
//============================================================================//
void class_superblock::find_ground_state(int eigSteps, double eigThreshold){
    Textra::Tensor1d theta = MPS.get_theta().shuffle(array4{1,0,2,3}).reshape(shape1);

    Spectra::SymEigsSolver<double,
            Spectra::SMALLEST_ALGE,
            class_superblock>
            eigs(this, 1, std::min(10,(int)theta.size()) );

    eigs.init(theta.data()); //Throw in the initial vector v0 here
    eigs.compute(eigSteps,eigThreshold, Spectra::SMALLEST_ALGE);

    if(eigs.info() != Spectra::SUCCESSFUL){
        cout << "Eigenvalue solver failed." << '\n';
        exit(1);
    }
    ground_state =  Matrix_to_Tensor<2,double>(eigs.eigenvectors(), shape2);

}

//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS.
//============================================================================//
void class_superblock::truncate(long chi_max_, double SVDThreshold){

    Eigen::JacobiSVD<MatrixXd> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(tensor2_to_matrix(ground_state), Eigen::ComputeThinU | Eigen::ComputeThinV);
    chi_max             = chi_max_;
    chi                 = std::min(SVD.rank(),chi_max);
//    double renorm       = (SVD.singularValues().head(chi).cwiseAbs2().sum());
    truncation_error    = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    U                   = Matrix_to_Tensor<3,double>(SVD.matrixU().leftCols(chi),{d,chiL,chi});
    MPS.LA              = Matrix_to_Tensor<1,double>(SVD.singularValues().head(chi).normalized(), {chi});
    V                   = Matrix_to_Tensor<3,double>(SVD.matrixV().leftCols(chi),{d,chiR,chi}).shuffle(array3{0,2,1});

}




//============================================================================//
// Canonicalization, or orthonormalization of the MPS
//============================================================================//


// According to Orus and Vidal
//void class_superblock::canonicalize_iMPS(){
//    if (MPS.L_tail.size() <= 1){
//        return;
//    }
//    Eigen::JacobiSVD<MatrixXd> SVD;
//    SelfAdjointEigenSolver<MatrixXd> esR;
//    SelfAdjointEigenSolver<MatrixXd> esL;
////    cout << "Swapped: " << MPS.swapped << '\n';
//    //Coarse grain
//    update_bond_dimensions(); //Gives us chiL __|__chi__|__chiR
//
//    // G is a Tensor4d  chia__|__|__chib where we expect chiL = chiR for infinite algorithms,
//    // but for finite ones chiL=/=chiR because GA.dimension(1) =/= GB.dimension(2).
//    Tensor3d GA     = MPS.GA;
//    Tensor3d GA_old = GA;
//    //L is simply of size chiR because L = LB;
//    Tensor1d LA = MPS.LA;
//    Tensor1d LB = MPS.LB;
//    long chiLA = LA.size();
//    long chiLB = LB.size();
//    Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});
//    Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1,0,2});
//
//    //Watch out for a spurious bug without .eval() here!
//    Tensor2d TR = GALA
//            .contract(GALA.conjugate(), idxlist1{idx2(0, 0)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiLB * chiLB, chiLA * chiLA});
//    Tensor2d TL = LBGA
//            .contract(LBGA.conjugate(), idxlist1{idx2(0, 0)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiLB * chiLB, chiLA * chiLA});
//    cout << setprecision(16);
//
//
//    cout << "LA old: \n" << MPS.LA << endl;
//    cout << "LB old: \n" << MPS.LB << endl;
////        return;
//    cout << "Transf R: \n" << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << endl;
//    cout << "Transf L: \n"<<  LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << endl;
//    cout << "TR\n" << TR << endl;
//    cout << "TL\n" << TL << endl;
//
//    // Wikipedia:
//    //Left and right eigenvectors[edit]
////    See also: left and right (algebra)
////    Many disciplines traditionally represent vectors as matrices with a single column rather than as matrices with a single row. For that reason, the word "eigenvector" in the context of matrices almost always refers to a right eigenvector, namely a column vector that right multiples the n by n matrix A in the defining equation, Equation (1),
////
////    {\displaystyle Av=\lambda v.} Av=\lambda v.
////            The eigenvalue and eigenvector problem can also be defined for row vectors that left multiply matrix A. In this formulation, the defining equation is
////
////            {\displaystyle uA=\kappa u,} {\displaystyle uA=\kappa u,}
////    where κ is a scalar and u is a 1 by n matrix. Any row vector u satisfying this equation is called a left eigenvector of A and κ is its associated eigenvalue. Taking the transpose of this equation,
////
////            {\displaystyle A^{T}u^{T}=\kappa u^{T}.} {\displaystyle A^{T}u^{T}=\kappa u^{T}.}
////    Comparing this equation to Equation (1), it follows immediately that a left eigenvector of A is the same as the transpose of a right eigenvector of AT, with the same eigenvalue. Furthermore, since the characteristic polynomial of AT is the same as the characteristic polynomial of A, the eigenvalues of the left eigenvectors of A are the same as the eigenvalues of the right eigenvectors of AT.
////
//    Tensor2d VR = Math::dominant_eigenvector<Math::Handedness::RGHT>(TR).real().reshape(array2{chiLB, chiLB});
//    Tensor2d VL = Math::dominant_eigenvector<Math::Handedness::RGHT >(TL.shuffle(array2{1,0})).real().reshape(array2{chiLA, chiLA});
////
////    esR.compute(Tensor2d_to_MatrixXd(TR));
////    esL.compute(Tensor2d_to_MatrixXd(TL));
////    Tensor2d VRe = MatrixXd_to_Tensor2d(esR.eigenvectors().rightCols(1)).reshape(array2{chiR, chiR});
////    Tensor2d VLe = MatrixXd_to_Tensor2d(esL.eigenvectors().rightCols(1)).reshape(array2{chiR, chiR});
//
//    cout << "VR    =  \n" << VR << '\n';
//    cout << "VL    =  \n" << VL << '\n';
//
//
//    //============================================================================//
//    // Find square root using SVD decomposition
//    //============================================================================//
////    cout << "\n\nSVD:\n";
////
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
//    cout << "\n\nEigenvalue:\n";
//    esR.compute(Tensor2d_to_MatrixXd(VR));
////    Tensor2d _X    = MatrixXd_to_Tensor2d( esR.eigenvectors()*esR.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt());
////    Tensor2d _XT   = MatrixXd_to_Tensor2d( esR.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt() * esR.eigenvectors().transpose());
////    Tensor2d _Xi   = Tensor2d_inverse(_X);
//    Tensor2d _X    = MatrixXd_to_Tensor2d(esR.operatorSqrt());
//    Tensor2d _Xi   = MatrixXd_to_Tensor2d(esR.operatorInverseSqrt());
////    cout << "X : \n" << _X << endl;
////    cout << "Xi : \n" << _Xi << endl;
//
////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//    esL.compute(Tensor2d_to_MatrixXd(VL));
//    Tensor2d _Y    = MatrixXd_to_Tensor2d(esL.operatorSqrt());
//    Tensor2d _Yi   = MatrixXd_to_Tensor2d(esL.operatorInverseSqrt());
////        cout << "YT : \n"  << _YT << endl;
////        cout << "YTi : \n" << _YTi << endl;
////        cout << "YTiYT: \n" << YTi.contract(YT, idxlist1{idx2(1,0)}) << endl;
////        cout << "YTY: \n"   << YT.contract(Y, idxlist1{idx2(1,0)}) << endl;
//
//
//    Tensor2d YLAX = _Y.contract(asDiagonal(LA), idxlist1{idx2(1,0)})
//            .contract(_X, idxlist1{idx2(1,0)});
//
//    Tensor2d YLBX = _Y.contract(asDiagonal(LB), idxlist1{idx2(1,0)})
//            .contract(_X, idxlist1{idx2(1,0)});
//
//
//
//    SVD.compute(Tensor2d_to_MatrixXd(YLAX), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU());
//    Tensor2d VA = MatrixXd_to_Tensor2d(SVD.matrixV().transpose());
////    cout << "U: \n" << U << endl;
////    cout << "V: \n" << V << endl;
//    LA         = MatrixXd_to_Tensor1d(SVD.singularValues().normalized());
//
//
//    SVD.compute(Tensor2d_to_MatrixXd(YLBX), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU());
//    Tensor2d VB = MatrixXd_to_Tensor2d(SVD.matrixV().transpose());
////    cout << "U: \n" << U << endl;
////    cout << "V: \n" << V << endl;
//    LB         = MatrixXd_to_Tensor1d(SVD.singularValues().normalized());
//
//
//
//    GA         = VA
//            .contract(_Xi, idxlist1{idx2(1,0)})
//            .contract(GA_old, idxlist1{idx2(1,1)})
//            .contract(_Yi, idxlist1{idx2(2, 0)})
//            .contract(UB, idxlist1{idx2(2,0)})
//            .shuffle(array3{1,0,2});
//
//
//    cout << "LA new: \n" << LA << endl;
//    cout << "LB new: \n" << LB << endl;
//
//    cout << "\n\n Sanity check \n" << endl;
//    GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});
//    LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1,0,2});
//
//    //Watch out for a spurious bug without .eval() here!
//    TR = GALA
//            .contract(GALA.conjugate(), idxlist1{idx2(0, 0)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiLB * chiLB, chiLA * chiLA});
//    TL = LBGA
//            .contract(LBGA.conjugate(), idxlist1{idx2(0, 0)})
//            .shuffle(array4{0, 2, 1, 3})
//            .reshape(array2{chiLB * chiLB, chiLA * chiLA});
//    cout << setprecision(16);
//
//
//    VR = Math::dominant_eigenvector<Math::Handedness::RGHT>(TR).real().reshape(array2{chiLB, chiLB});
//    VL = Math::dominant_eigenvector<Math::Handedness::LEFT>(TL).real().reshape(array2{chiLA, chiLA});
//
//    cout << "VR    =  \n" << VR << '\n';
//    cout << "VL    =  \n" << VL << '\n';
//    cout << "TR    =  \n" << TR << '\n';
//    cout << "TL    =  \n" << TL << '\n';
//    cout << "Transf R: \n" << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << endl;
//    cout << "Transf L: \n"<<  LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << endl;
//
//    cout << "Eigen diagonalizer:\n";
//    esR.compute(Tensor2d_to_MatrixXd(TR));
//    esL.compute(Tensor2d_to_MatrixXd(TL));
//    cout << "VR eigval   =  \n" << esR.eigenvalues() << '\n';
//    cout << "VR eigvec   =  \n" << esR.eigenvectors() << '\n';
//    cout << "VL eigval   =  \n" << esL.eigenvalues() << '\n';
//    cout << "VL eigvec   =  \n" << esL.eigenvectors() << '\n';
//
//
//    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib
//
//}



// According to Orus and Vidal
void class_superblock::canonicalize2(){
    if (MPS.L_tail.size() <= 1){
        return;
    }
    Eigen::JacobiSVD<MatrixXd> SVD;
    SelfAdjointEigenSolver<MatrixXd> esR;
    SelfAdjointEigenSolver<MatrixXd> esL;
    //Coarse grain
    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib

    // G is a Tensor4d  chia__|__|__chib where we expect chiL = chiR for infinite algorithms,
    // but for finite ones chiL=/=chiR because GA.dimension(1) =/= GB.dimension(2).
    Tensor4d G = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)})
            .contract(MPS.GB, idxlist1{idx2(2, 1)}).shuffle(array4{0, 2, 1, 3});
    Tensor4d G_old = G;
    //L is simply of size chiR because L = LB;
    Tensor1d L = MPS.LB;

    Tensor4d GL = G.contract(asDiagonal(L), idxlist1{idx2(3, 0)});
    Tensor4d LG = asDiagonal(L).contract(G, idxlist1{idx2(1, 2)}).shuffle(array4{1,2,0,3});

    //Watch out for a spurious bug without .eval() here!
    Tensor2d TR = GL
            .contract(GL.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
            .shuffle(array4{0, 2, 1, 3})
            .reshape(array2{chiR * chiR, chiR * chiR});
    Tensor2d TL = LG
            .contract(LG.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
            .shuffle(array4{0, 2, 1, 3})
            .reshape(array2{chiR * chiR, chiR * chiR});
    cout << setprecision(16);


    cout << "LA old: \n" << MPS.LA << endl;
    cout << "LB old: \n" << MPS.LB << endl;
    cout << "TR: \n" << GL.contract(GL.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(3,3)}) << endl;
    cout << "TL: \n"<<  LG.contract(LG.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(2,2)}) << endl;

    Tensor2d VR = Math::dominant_eigenvector<Math::Handedness::RGHT>(TR).real().reshape(array2{chiR, chiR});
    Tensor2d VL = Math::dominant_eigenvector<Math::Handedness::RGHT>(TL.shuffle(array2{1,0})).real().reshape(array2{chiR, chiR});

    cout << "VR    =  \n" << VR << '\n';
    cout << "VL    =  \n" << VL << '\n';


    //============================================================================//
    // Find square root using SVD decomposition
    //============================================================================//
//    SVD.compute(Tensor2d_to_MatrixXd(VR), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    Tensor2d X               = MatrixXd_to_Tensor2d( SVD.matrixU()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixU().transpose());
//    Tensor2d Xi              = Tensor2d_inverse(X);
//    cout << "X : \n" << X << endl;
//    cout << "Xi : \n" << Xi << endl;
////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//    SVD.compute(Tensor2d_to_MatrixXd(VL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    Tensor2d Y               = MatrixXd_to_Tensor2d( SVD.matrixV().transpose()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixV());
//    Tensor2d Yi              = Tensor2d_inverse(Y);
//    cout << "Y : \n" << Y << endl;
//    cout << "Yi : \n" << Yi << endl;
//        cout << "YTiYT: \n" << YTi.contract(YT, idxlist1{idx2(1,0)}) << endl;
//        cout << "YTY: \n" << YT.contract(Y, idxlist1{idx2(1,0)}) << endl;


//    ============================================================================//
//     Find square root using Eigenvalue decomposition
//    ============================================================================//
    esR.compute(Tensor2d_to_MatrixXd(VR));
    Tensor2d _X    = MatrixXd_to_Tensor2d(esR.operatorSqrt());
    Tensor2d _XT   = _X.shuffle(array2{1,0});
    Tensor2d _Xi   = MatrixXd_to_Tensor2d(esR.operatorInverseSqrt());
    cout << "XXT: \n" << _X.contract(_XT, idxlist1{idx2(1,0)}) << endl;

    esL.compute(Tensor2d_to_MatrixXd(VL));
    Tensor2d _YT    = MatrixXd_to_Tensor2d(esL.operatorSqrt().transpose());
    Tensor2d _Y     = MatrixXd_to_Tensor2d(esL.operatorSqrt());
    Tensor2d _YTi   = MatrixXd_to_Tensor2d(esL.operatorInverseSqrt().transpose());
    cout << "YTY: \n"   << _YT.contract(_Y, idxlist1{idx2(1,0)}) << endl;


    Tensor2d YLX  = _YT.contract(asDiagonal(L), idxlist1{idx2(1,0)})
            .contract(_X, idxlist1{idx2(1,0)});



    SVD.compute(Tensor2d_to_MatrixXd(YLX), Eigen::ComputeThinU | Eigen::ComputeThinV);
    Tensor2d U = MatrixXd_to_Tensor2d(SVD.matrixU());
    Tensor2d V = MatrixXd_to_Tensor2d(SVD.matrixV().transpose());
    L          = MatrixXd_to_Tensor1d(SVD.singularValues().normalized());

    G         = V
            .contract(_Xi, idxlist1{idx2(1,0)})
            .contract(G_old, idxlist1{idx2(1,2)})
            .contract(_YTi, idxlist1{idx2(3, 0)})
            .contract(U, idxlist1{idx2(3,0)})
            .shuffle(array4{1,2,0,3});

    cout << "L new: \n" << L << endl;


    cout << "\n\n Sanity check \n" << endl;
    GL = G.contract(asDiagonal(L), idxlist1{idx2(3, 0)});
    LG = G.contract(asDiagonal(L), idxlist1{idx2(2, 1)}).shuffle(array4{0, 1, 3, 2});
    //Watch out for a spurious bug without .eval() here!
    TR = GL
            .contract(GL.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
            .shuffle(array4{0, 2, 1, 3})
            .reshape(array2{chiR * chiR, chiR * chiR});
    TL = LG
            .contract(LG.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)})
            .shuffle(array4{0, 2, 1, 3})
            .reshape(array2{chiR * chiR, chiR * chiR});
    cout << "TR    =  \n" << TR << '\n';
    cout << "TL    =  \n" << TL << '\n';
    cout << "Transf R: \n" << GL.contract(GL.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(3,3)}) << endl;
    cout << "Transf L: \n"<<  LG.contract(LG.conjugate(), idxlist3{idx2(0, 0), idx2(1, 1), idx2(2,2)}) << endl;
    VR = Math::dominant_eigenvector<Math::Handedness::RGHT>(TR).real().reshape(array2{chiR, chiR});
    VL = Math::dominant_eigenvector<Math::Handedness::RGHT>(TL.shuffle(array2{1,0})).real().reshape(array2{chiR, chiR});

    cout << "VR    =  \n" << VR << '\n';
    cout << "VL    =  \n" << VL << '\n';






    //Separate
    Tensor2d LGL =   asDiagonal(L)
            .contract(G, idxlist1{idx2(1,2)})
            .contract(asDiagonal(L), idxlist1{idx2(3,0)})
            .shuffle(array4{1,0,2,3})
            .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});

    SVD.compute(Tensor2d_to_MatrixXd(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
    chiL                         = std::min(SVD.rank(),chi_max);
    truncation_error             = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiL).sum();
    Tensor3d U3                   = Matrix_to_Tensor<3,double>(SVD.matrixU(),{d,chiR,chiL});
    MPS.LA                        = Matrix_to_Tensor<1,double>(SVD.singularValues().head(chiL).normalized(), {chiL});
    Tensor3d V3                   = Matrix_to_Tensor<3,double>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});

    MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
    MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
    MPS.LB      = Tensor1d_normalize(L);
    MPS.L_tail  = MPS.LB;
    cout << "LA new: \n" << MPS.LA << endl;
    cout << "LB new: \n" << MPS.LB << endl;
//        MPS.swap_AB();


    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib

}



// According to Ran, Shi-Ju, Emanuele Tirrito, Cheng Peng, Xi Chen, Gang Su, and Maciej Lewenstein. 2017. “Review of Tensor Network Contraction Approaches,” September. http://arxiv.org/abs/1708.09213.
void class_superblock::canonicalize4(){
    if (MPS.L_tail.size() <= 1){
        return;
    }
    cout << setprecision(16);

    cout << "LA old: \n" << MPS.LA << endl;
    cout << "LB old: \n" << MPS.LB << endl;
    SelfAdjointEigenSolver<MatrixXd> esR;
    SelfAdjointEigenSolver<MatrixXd> esL;
    Eigen::JacobiSVD<MatrixXd> SVD;
    SVD.setThreshold(settings::precision::SVDThreshold);
    for(int i = 0; i < 2; i++ ) {
        update_bond_dimensions(); //Gives us chiL __|__chi__|__chiR
        //Following the article:
        //Ran, Shi-Ju, Emanuele Tirrito, Cheng Peng, Xi Chen, Gang Su, and Maciej Lewenstein. 2017. “Review of Tensor Network Contraction Approaches,” September. http://arxiv.org/abs/1708.09213.
        // Lambda -> Lambda^B
        // Gamma  -> Lambda^A
        // A      -> Gamma^A
        // B      -> Gamma^B
        long chiA = MPS.LA.size();
        long chiB = MPS.LB.size();
        Tensor3d LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
        Tensor3d GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
        Tensor3d LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
        Tensor3d GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
        Tensor4d ML = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
        Tensor4d MR = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
        Tensor4d NL = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
        Tensor4d NR = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});

        cout << "ML dim: " << ML.dimensions() << endl;
        cout << "MR dim: " << MR.dimensions() << endl;
        cout << "NL dim: " << NL.dimensions() << endl;
        cout << "NR dim: " << NR.dimensions() << endl;
//
//        Tensor2d NRML = NR.contract(ML, idxlist2{idx2(2, 0), idx2(3, 1)})
//                           .shuffle(array4{0,1,2,3})
//                           .reshape(array2{chiA*chiA,chiA*chiA}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)
//        Tensor2d NLMR = NL.contract(MR, idxlist2{idx2(2, 0), idx2(3, 1)})
//                          .shuffle(array4{0,1,2,3})
//                          .reshape(array2{chiA*chiA,chiA*chiA}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)

        Tensor2d NRML = ML.contract(NR, idxlist2{idx2(2, 0), idx2(3, 1)})
                .shuffle(array4{0,1,2,3})
                .reshape(array2{chiB*chiB,chiB*chiB}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)
        Tensor2d NLMR = MR.contract(NL, idxlist2{idx2(2, 0), idx2(3, 1)})
                .shuffle(array4{0,1,2,3})
                .reshape(array2{chiB*chiB,chiB*chiB}); // (chi*chi , chiR*chiR) X (chiL*chiL , chi*chi)


//        cout << "NRML: \n" << NRML << '\n';
//        cout << "NLMR: \n" << NLMR << '\n';
        Tensor2d VR = Math::dominant_eigenvector<Math::Handedness::RGHT>  (NRML).real().reshape(array2{chiB, chiB});
        Tensor2d VL = Math::dominant_eigenvector<Math::Handedness::RGHT > (NLMR.shuffle(array2{1,0})).real().reshape(array2{chiB, chiB});
        cout << "VR: \n" << VR << '\n';
        cout << "VL: \n" << VL << '\n';

        //    ============================================================================//
        //     Find square root using Eigenvalue decomposition
        //    ============================================================================//
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


        SVD.compute(Tensor2d_to_MatrixXd(VR), Eigen::ComputeThinU | Eigen::ComputeThinV);
        Tensor2d X = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixV().transpose());
        Tensor2d XT = X.shuffle(array2{1, 0});
        Tensor2d Xi = Tensor2d_inverse(X);
        SVD.compute(Tensor2d_to_MatrixXd(VL), Eigen::ComputeThinU | Eigen::ComputeThinV);
        Tensor2d Y = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixV().transpose());
        Tensor2d YT = Y.shuffle(array2{1, 0});
        Tensor2d YTi = Tensor2d_inverse(YT);

        cout << "X      =  \n" << X << '\n';
        cout << "X*XT   =  \n" << X.contract(XT, idxlist1{idx2(1, 0)}) << '\n';
        cout << "Xi*X   =  \n" << Xi.contract(X, idxlist1{idx2(1, 0)}) << '\n';
        cout << "Y      =  \n" << Y << '\n';
        cout << "YT*Y   =  \n" << YT.contract(Y, idxlist1{idx2(1, 0)}) << '\n';
        cout << "YT*YTi   =  \n" << YTi.contract(YT, idxlist1{idx2(1, 0)}) << '\n';


        Tensor2d X_LA_Y = X
                .contract(asDiagonal(MPS.LA), idxlist1{idx2(1, 0)})
                .contract(Y, idxlist1{idx2(1, 0)});

        SVD.compute(Tensor2d_to_MatrixXd(X_LA_Y), Eigen::ComputeThinU | Eigen::ComputeThinV);
        chiA = std::min(SVD.rank(),chi_max);
        Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
        Tensor2d VA = MatrixXd_to_Tensor2d(SVD.matrixV().leftCols(chiA).transpose());
        MPS.LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());


        Tensor3d GA = MPS.GA.contract(Xi, idxlist1{idx2(2, 0)})
                .contract(UA, idxlist1{idx2(2, 0)});
//        Tensor3d GB = VA.contract(YTi, idxlist1{idx2(1, 0)}) //Transpose Y?
//                .contract(MPS.GB, idxlist1{idx2(1, 1)})
//                .shuffle(array3{1, 0, 2});
        Tensor3d GB = MPS.GB.contract(YTi, idxlist1{idx2(1, 1)})
                .contract(VA, idxlist1{idx2(2, 1)})
                .shuffle(array3{0,2,1});
        MPS.GA = GA;
        MPS.GB = GB;
        cout << "LA new: \n" << MPS.LA << endl;
        cout << "LB new: \n" << MPS.LB << endl;
        MPS.L_tail = MPS.LB;


        cout << "\n\n Sanity Check " << endl;
        LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
        GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
        LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
        GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
        Tensor2d ML_ = LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1,1)});
        Tensor2d MR_ = GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2,2)});
        Tensor2d NL_ = LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1,1)});
        Tensor2d NR_ = GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2,2)});
        cout << "ML: \n " <<ML_<< '\n';
        cout << "MR: \n " <<MR_<< '\n';
        cout << "NL: \n " <<NL_<< '\n';
        cout << "MR: \n " <<MR_<< '\n';






        MPS.swap_AB();
    }



    update_bond_dimensions(); //Gives us chiL __|__chi__|__chib
}


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
//        SVD.compute(tensor2_to_matrix<double>(VR), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        Tensor2d X               = matrix_to_tensor2<double>( SVD.matrixU()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixU().transpose());
//        Tensor2d XT              = X.shuffle(array2{1,0});
//        Tensor2d Xi              = matrix_to_tensor2<double>(tensor2_to_matrix(X).inverse());
//        cout << "X : \n" << X << endl;
//        cout << "Xi : \n" << Xi << endl;
////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//        SVD.compute(tensor2_to_matrix(VL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//        Tensor2d Y               = matrix_to_tensor2<double>( SVD.matrixU()*SVD.singularValues().cwiseSqrt().asDiagonal()*SVD.matrixU().transpose());
//        Tensor2d Yi              = matrix_to_tensor2<double>(tensor2_to_matrix(Y).inverse());
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
//        esR.compute(tensor2_to_matrix(VR));
//        Tensor2d _X    = matrix_to_tensor2<double>( esR.eigenvectors()*esR.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt());
//        Tensor2d _XT   = matrix_to_tensor2<double>( esR.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt() * esR.eigenvectors().transpose());
//        Tensor2d _Xi   = matrix_to_tensor2<double>(tensor2_to_matrix(X).inverse());
//        cout << "X : \n" << _X << endl;
//        cout << "Xi : \n" << _Xi << endl;
////        cout << "XXi: \n" << X.contract(Xi, idxlist1{idx2(1,0)}) << endl;
////        cout << "XXT: \n" << X.contract(XT, idxlist1{idx2(1,0)}) << endl;
//
//        esL.compute(tensor2_to_matrix(VL));
//        Tensor2d _Y              = matrix_to_tensor2<double>( esL.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt() * esL.eigenvectors().transpose());
//        Tensor2d _YT             = matrix_to_tensor2<double>( esL.eigenvectors()*esL.eigenvalues().asDiagonal().toDenseMatrix().cwiseSqrt());
//        Tensor2d _YTi            = matrix_to_tensor2<double>(tensor2_to_matrix(_YT).inverse());
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
//        SVD.compute(tensor2_to_matrix(YLX), Eigen::ComputeThinU | Eigen::ComputeThinV);
//
//        Tensor2d U = matrix_to_tensor2<double>(SVD.matrixU());
//        Tensor2d V = matrix_to_tensor2<double>(SVD.matrixV().transpose());
//        L         = Matrix_to_Tensor1<double>(SVD.singularValues().normalized());
//
//        G         = V
//                   .contract(Xi, idxlist1{idx2(1,0)})
//                   .contract(G_old, idxlist1{idx2(1,2)})
//                   .contract(Yi, idxlist1{idx2(3, 0)})
//                   .contract(U, idxlist1{idx2(3,0)})
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
//        SVD.compute(tensor2_to_matrix(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
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

