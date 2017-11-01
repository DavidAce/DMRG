//
// Created by david on 2017-10-05.
//
#include <class_superblock.h>
#include <n_math.h>
#include <iomanip>
#include <SymEigsSolver.h>
#include <Eigen/Eigenvalues>


void class_superblock::test(){
    auto tuple = std::make_tuple(1, 'a', 2.3);

    // unpack the tuple into individual variables declared at the call site
    auto [ i, c, d ] = tuple;

    std::cout << "i=" << i << " c=" << c << " d=" << d << '\n';
    return;
}

void class_superblock::transfer_matrix_eigenvalues() {

    if (MPS.LB.size()!= MPS.LA.size()){
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    if (MPS.LB.size() < 2 || MPS.LA.size() < 2 ){
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    cout << setprecision(16);
    update_bond_dimensions();
    long chiA = MPS.LA.size();
    long chiB = MPS.LB.size();
    Tensor3d GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
    Tensor3d GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
    Tensor3d LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
    Tensor3d LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)

    Tensor4d TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});

    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
    Tensor2d VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
    Tensor2d VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
    Spectra::DenseGenMatProd<double> op_TA_R(Tensor_to_MatrixXd<4>(TA_R, chiB*chiB,chiA*chiA));
    Spectra::DenseGenMatProd<double> op_TA_L(Tensor_to_MatrixXd<4>(TA_L, chiB*chiB,chiA*chiA));
    Spectra::DenseGenMatProd<double> op_TB_R(Tensor_to_MatrixXd<4>(TB_R, chiA*chiA,chiB*chiB));
    Spectra::DenseGenMatProd<double> op_TB_L(Tensor_to_MatrixXd<4>(TB_L, chiA*chiA,chiB*chiB));
    int ncv = std::min(settings::precision::eig_max_ncv,4);
    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTA_R(&op_TA_R, 2, ncv);
    eigsTA_R.init();
    eigsTA_R.compute(10000, 1e-12, Spectra::LARGEST_REAL);
    std::cout << "Eigenvalue TA_R:\n" << eigsTA_R.eigenvalues()<< std::endl;
    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTA_L(&op_TA_L, 2, ncv);
    eigsTA_L.init();
    eigsTA_L.compute(10000, 1e-12, Spectra::LARGEST_REAL);
    std::cout << "Eigenvalue TA_L:\n" << eigsTA_L.eigenvalues()<< std::endl;
    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTB_R(&op_TB_R, 2, ncv);
    eigsTB_R.init();
    eigsTB_R.compute(10000, 1e-12, Spectra::LARGEST_REAL);
    std::cout << "Eigenvalue TB_R:\n" << eigsTB_R.eigenvalues()<< std::endl;
    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>> eigsTB_L(&op_TB_L, 2, ncv);
    eigsTB_L.init();
    eigsTB_L.compute(10000, 1e-12, Spectra::LARGEST_REAL);
    std::cout << "Eigenvalue TB_L: \n" << eigsTB_L.eigenvalues()<< std::endl;

    return;
}

void class_superblock::canonicalize_iMPS_iterative_again() {
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
    cout << setprecision(16) << endl;
    cout << "LA old: \n" << MPS.LA << endl;
    cout << "LB old: \n" << MPS.LB << endl;
    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;

    Eigen::BDCSVD<MatrixXd> SVD;
//    SelfAdjointEigenSolver<MatrixXd> es;
//
    SVD.setThreshold(1e-24);
    int iter = 0;
    Tensor3d GA = MPS.GA;
    Tensor3d GB = MPS.GB;
    Tensor1d LA = MPS.LA;
    Tensor1d LB = MPS.LB;
    update_bond_dimensions();

    while (diffA + diffB > 1e-10) {
//        cout << "diff: " << diff << endl;
        iter++;
        //Coarse grain
        Tensor3d GA_old = GA;
        Tensor3d GB_old = GB;
        Tensor1d LA_old = LA;
        Tensor1d LB_old = LB;

        long chiA = LA.size();
        long chiB = LB.size();
        Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
        Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
        Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
        Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)

        Tensor2d VA_R = GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)});
        Tensor2d VA_L = LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)});
        Tensor2d VB_R = GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)});
        Tensor2d VB_L = LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)});

        //============================================================================//
        // Find square root using Eigenvalue decomposition
        //============================================================================//
        SelfAdjointEigenSolver<MatrixXd> es;
        es.compute(Tensor2d_to_MatrixXd(VA_R));
        Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
        Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
        es.compute(Tensor2d_to_MatrixXd(VA_L));
        Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
        Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
        es.compute(Tensor2d_to_MatrixXd(VB_R));
        Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
        Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
        es.compute(Tensor2d_to_MatrixXd(VB_L));
        Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
        Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());

        //Symmetric
        Tensor2d YA_LA_XB = YA
                .contract(asDiagonal(LA), idxlist1{idx2(1, 0)})
                .contract(XB, idxlist1{idx2(1, 0)});
        //Symmetric
        Tensor2d YB_LB_XA = YB
                .contract(asDiagonal(LB), idxlist1{idx2(1, 0)})
                .contract(XA, idxlist1{idx2(1, 0)});

        SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
        chiA = std::min(SVD.rank(), chi_max);
        Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
        Tensor2d VA = MatrixXd_to_Tensor2d(
                SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, U and V are not symmetric!
        LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());

        SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
        chiB = std::min(SVD.rank(), chi_max);
        Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
        Tensor2d VB = MatrixXd_to_Tensor2d(
                SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, U and V are not symmetric!
        LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());

//    cout << "UA: \n "<< UA << endl;
//    cout << "VA: \n "<< VA << endl;
//    cout << "UB: \n "<< UB << endl;
//    cout << "VB: \n "<< VB << endl;
        GA = VB
                .contract(XAi, idxlist1{idx2(1, 0)})
                .contract(GA_old, idxlist1{idx2(1, 1)})
                .contract(YAi, idxlist1{idx2(2, 0)})
                .contract(UA, idxlist1{idx2(2, 0)})
                .shuffle(array3{1, 0, 2});
        GB = VA
                .contract(XBi, idxlist1{idx2(1, 0)})
                .contract(GB_old, idxlist1{idx2(1, 1)})
                .contract(YBi, idxlist1{idx2(2, 0)})
                .contract(UB, idxlist1{idx2(2, 0)})
                .shuffle(array3{1, 0, 2});


        Tensor0d tempA = (LA - LA_old.slice(array1{0}, array1{chiA})).square().sum().sqrt();
        Tensor0d tempB = (LB - LB_old.slice(array1{0}, array1{chiB})).square().sum().sqrt();
        diffA = tempA(0);
        diffB = tempB(0);
//         cout << "LA diff: " << diffA << endl;
//         cout << "LB diff: " << diffB << endl;
//        if(iter > 100 ){
//            exit(0);
//        }
    }
    Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
    Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
    Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
    Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
    cout << "TA_R: \n " << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
    cout << "TA_L: \n " << LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
    cout << "TB_R: \n " << GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
    cout << "TB_L: \n " << LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
//    Tensor2d LGL =   asDiagonal(L)
//                   .contract(G, idxlist1{idx2(1,2)})
//                   .contract(asDiagonal(L), idxlist1{idx2(3,0)})
//                   .shuffle(array4{1,0,2,3})
//                   .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});
//
//
////    //Note the difference! U3 is a d x [chib] x chi matrix! Previously we had [chia] in normal SVD truncation
//    SVD.compute(tensor2_to_matrix(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    chiL                         = std::min(SVD.rank(),chi_max);
//    truncation_error             = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiL).sum();
//    Tensor3d U3                   = matrix_to_tensor<3>(SVD.matrixU(),{d,chiR,chiL});
//    MPS.LA                       = matrix_to_tensor<1>(SVD.singularValues().head(chiL).normalized(), {chiL});
//    Tensor3d V3                   = matrix_to_tensor<3>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});
////
//    MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
//    MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
//    MPS.LB = L;
//    MPS.L_tail = L;

//    cout << "Full LA after: \n" << SVD.singularValues().normalized() << endl;
    cout << "LA after: \n" << LA << endl;
    cout << "LB after: \n" << LB << endl;
    cout << "Entropy A: " << -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
    cout << "Entropy B: " << -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
    cout << "\n\niter: " << iter << endl;

//     MPS.GA = GA;
//     MPS.GB = GB;
//     MPS.LA = LA;
//     MPS.LB = LB;
//     MPS.L_tail = LB;
//    TRY USING ROBUST SVD!
//             Stoudenmire, E. M., & White, S. R. (2013). Real-space parallel density matrix renormalization group. Physical Review B, 87(April), 155137. https://doi.org/10.1103/PhysRevB.87.155137

}

void class_superblock::canonicalize_iMPS_vidal() {
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
//    double diffA = 1;
//    double diffB = 1;
    cout << setprecision(16) << endl;
    cout << "LA old: \n" << MPS.LA << endl;
    cout << "LB old: \n" << MPS.LB << endl;
    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;

    Eigen::BDCSVD<MatrixXd> SVD;
//    SelfAdjointEigenSolver<MatrixXd> es;
//
    SVD.setThreshold(1e-24);
    int iter = 0;
    Tensor3d GA = MPS.GA;
    Tensor3d GB = MPS.GB;
    Tensor1d LA = MPS.LA;
    Tensor1d LB = MPS.LB;
    update_bond_dimensions();

    Tensor3d GA_old = GA;
    Tensor3d GB_old = GB;
    Tensor1d LA_old = LA;
    Tensor1d LB_old = LB;

    long chiA = LA.size();
    long chiB = LB.size();
    Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
    Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
    Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
    Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)

    Tensor4d TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});


    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
    Tensor2d VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
    Tensor2d VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});

    VA_R = VA_R/VA_R.coeff(0);
    VA_L = VA_L/VA_L.coeff(0);
    VB_R = VB_R/VB_R.coeff(0);
    VB_L = VB_L/VB_L.coeff(0);
    cout << "VA_R: \n " << VA_R << '\n';
    cout << "VA_L: \n " << VA_L << '\n';
    cout << "VB_R: \n " << VB_R << '\n';
    cout << "VB_L: \n " << VB_L << '\n';
    //============================================================================//
    // Find square root using Eigenvalue decomposition
    //============================================================================//
    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(Tensor2d_to_MatrixXd(VA_R));
    Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VA_L));
    Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VB_R));
    Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VB_L));
    Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());

    //Symmetric
    Tensor2d YA_LA_XB = YA
            .contract(asDiagonal(LA), idxlist1{idx2(1, 0)})
            .contract(XB, idxlist1{idx2(1, 0)});
    //Symmetric
    Tensor2d YB_LB_XA = YB
            .contract(asDiagonal(LB), idxlist1{idx2(1, 0)})
            .contract(XA, idxlist1{idx2(1, 0)});

    SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
    chiA = std::min(SVD.rank(), chi_max);
    Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
    Tensor2d VA = MatrixXd_to_Tensor2d(
            SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, U and V are not symmetric!
    LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());

    SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
    chiB = std::min(SVD.rank(), chi_max);
    Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
    Tensor2d VB = MatrixXd_to_Tensor2d(
            SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, U and V are not symmetric!
    LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());

//    cout << "UA: \n "<< UA << endl;
//    cout << "VA: \n "<< VA << endl;
//    cout << "UB: \n "<< UB << endl;
//    cout << "VB: \n "<< VB << endl;
    GA = VB
            .contract(XAi, idxlist1{idx2(1, 0)})
            .contract(GA_old, idxlist1{idx2(1, 1)})
            .contract(YAi, idxlist1{idx2(2, 0)})
            .contract(UA, idxlist1{idx2(2, 0)})
            .shuffle(array3{1, 0, 2});
    GB = VA
            .contract(XBi, idxlist1{idx2(1, 0)})
            .contract(GB_old, idxlist1{idx2(1, 1)})
            .contract(YBi, idxlist1{idx2(2, 0)})
            .contract(UB, idxlist1{idx2(2, 0)})
            .shuffle(array3{1, 0, 2});


    Tensor0d tempA = (LA - LA_old.slice(array1{0}, array1{chiA})).square().sum().sqrt();
    Tensor0d tempB = (LB - LB_old.slice(array1{0}, array1{chiB})).square().sum().sqrt();

//         cout << "LA diff: " << diffA << endl;
//         cout << "LB diff: " << diffB << endl;
//        if(iter > 100 ){
//            exit(0);
//        }


    cout << "LA after: \n" << LA << endl;
    cout << "LB after: \n" << LB << endl;
    cout << "Entropy A: " << -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
    cout << "Entropy B: " << -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;


    GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
    LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
    GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
    LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)

    TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});

    VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
    VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
    VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
    VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});
    cout << "VA_R: \n " << VA_R << '\n';
    cout << "VA_L: \n " << VA_L << '\n';
    cout << "VB_R: \n " << VB_R << '\n';
    cout << "VB_L: \n " << VB_L << '\n';

    cout << "TA_R: \n " << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
    cout << "TA_L: \n " << LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
    cout << "TB_R: \n " << GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
    cout << "TB_L: \n " << LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
    cout << "\n\niter: " << iter << endl;

     MPS.GA = GA;
     MPS.GB = GB;
     MPS.LA = LA;
     MPS.LB = LB;
     MPS.L_tail = LB;
//    TRY USING ROBUST SVD!
//             Stoudenmire, E. M., & White, S. R. (2013). Real-space parallel density matrix renormalization group. Physical Review B, 87(April), 155137. https://doi.org/10.1103/PhysRevB.87.155137

}


void class_superblock::canonicalize_iMPS_mcculloch() {
    //Start by creating a 2-site overlap

    if (MPS.LB.size()!= MPS.LA.size()){
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    if (MPS.LB.size() < 2 || MPS.LA.size() < 2 ){
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    Eigen::BDCSVD<MatrixXd> SVD;
    SVD.setThreshold(1e-14);
    cout << setprecision(16);
    update_bond_dimensions();
    long chiA = MPS.LA.size();
    long chiB = MPS.LB.size();
    cout << "LA old: \n" << MPS.LA << endl;
    cout << "LB old: \n" << MPS.LB << endl;
    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;

    Tensor2d P      = asDiagonal(MPS.LB).contract(asInverseDiagonal(MPS.LB), idxlist1{idx2(1,0)});
    Tensor2d Q      = asInverseDiagonal(MPS.LB).contract(asDiagonal(MPS.LB), idxlist1{idx2(1,0)});
    Tensor4d thetaP =  asDiagonal(MPS.LB)
                      .contract(MPS.GA, idxlist1{idx2(1,1)})
                      .contract(asDiagonal(MPS.LA), idxlist1{idx2(2,0)})
                      .contract(MPS.GB, idxlist1{idx2(2,1)})
                      .contract(P, idxlist1{idx2(3,0)});
    Tensor4d Qtheta =  Q
            .contract(MPS.GA, idxlist1{idx2(1,1)})
            .contract(asDiagonal(MPS.LA), idxlist1{idx2(2,0)})
            .contract(MPS.GB, idxlist1{idx2(2,1)})
            .contract(asDiagonal(MPS.LB), idxlist1{idx2(3,0)});

    Tensor4d TL     = thetaP.contract(thetaP,idxlist2{idx2(1,1),idx2(2,2)}).shuffle(array4{0,2,1,3});
    Tensor4d TR     = Qtheta.contract(Qtheta,idxlist2{idx2(1,1),idx2(2,2)}).shuffle(array4{0,2,1,3});

    Tensor2d VL     = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TL.reshape(array2{chiB*chiB,chiB*chiB})).real().reshape(array2{chiB, chiB});
    Tensor2d VR     = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TR.reshape(array2{chiA*chiA,chiA*chiA})).real().reshape(array2{chiA, chiA});

    //Rescale

    cout << "VL : \n" << VL << endl;
    cout << "VR : \n" << VR << endl;
    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(Tensor2d_to_MatrixXd(VL));
    Tensor2d X  = MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d Xi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VR));
    Tensor2d Y  = MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d Yi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    cout << "X : \n" << X << endl;
    cout << "Y : \n" << Y << endl;

    Tensor2d EL = Xi.contract(VL,idxlist1{idx2(1,0)}).contract(Xi,idxlist1{idx2(1,0)});
    cout << "EL: \n" << EL << endl;

    Tensor2d LB = X.contract(asDiagonal(MPS.LB), idxlist1{idx2(1,0)}).contract(Y,idxlist1{idx2(1,0)});
    cout << "LB: \n" << LB << endl;
    cout << "S_E: " <<  -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;


}


void class_superblock::canonicalize_iMPS() {
    if (MPS.LB.size()!= MPS.LA.size()){
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    if (MPS.LB.size() < 2 || MPS.LA.size() < 2 ){
        //The left/right eigenvectors are only hermitian if LA and LB have the same dimension!
        return;
    }
    cout << setprecision(16);
    update_bond_dimensions();
    long chiA = MPS.LA.size();
    long chiB = MPS.LB.size();
    cout << "LA old: \n" << MPS.LA << endl;
    cout << "LB old: \n" << MPS.LB << endl;
    cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
    cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;

    Eigen::BDCSVD<MatrixXd> SVD;
    SVD.setThreshold(1e-14);
    Tensor3d GALA = MPS.GA.contract(asDiagonal(MPS.LA), idxlist1{idx2(2, 0)});                          //(eq 107)
    Tensor3d GBLB = MPS.GB.contract(asDiagonal(MPS.LB), idxlist1{idx2(2, 0)});                          //(eq 109)
    Tensor3d LBGA = asDiagonal(MPS.LB).contract(MPS.GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
    Tensor3d LAGB = asDiagonal(MPS.LA).contract(MPS.GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)

    Tensor4d TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    Tensor4d TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
//
//    Tensor4d MRNR = TA_R.contract(TB_R.conjugate(), idxlist2{idx2(2, 0), idx2(3,1)});
//    Tensor4d MLNL = TA_L.contract(TB_L.conjugate(), idxlist2{idx2(2, 2), idx2(1,3)});
////    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RIGHT>  (MRNR.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
////    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::RIGHT> (MLNL.reshape(array2{chiB*chiB,chiA*chiA}).shuffle(array2{1,0})).real().reshape(array2{chiA, chiA});
////

    Tensor2d VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
    Tensor2d VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
    Tensor2d VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
    Tensor2d VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});

    cout << "VA_R : \n" << VA_R << endl;
    cout << "VA_L : \n" << VA_L << endl;
    cout << "VB_R : \n" << VB_R << endl;
    cout << "VB_L : \n" << VB_L << endl;

//    SVD.compute(Tensor2d_to_MatrixXd(VA_R), Eigen::ComputeThinU);
//    Tensor2d XA  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    Tensor2d XAi = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    SVD.compute(Tensor2d_to_MatrixXd(VA_L), Eigen::ComputeThinU);
//    Tensor2d YA   = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    Tensor2d YAi  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    SVD.compute(Tensor2d_to_MatrixXd(VB_R), Eigen::ComputeThinU);
//    Tensor2d XB   = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    Tensor2d XBi  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    SVD.compute(Tensor2d_to_MatrixXd(VB_L), Eigen::ComputeThinU);
//    Tensor2d YB   = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//    Tensor2d YBi  = MatrixXd_to_Tensor2d(SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//

    //============================================================================//
    // Find square root using Eigenvalue decomposition
    //============================================================================//
    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(Tensor2d_to_MatrixXd(VA_R));
    Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VA_L));
    Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VB_R));
    Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
    es.compute(Tensor2d_to_MatrixXd(VB_L));
    Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
    Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());

//    cout << "XA         =  \n"  << XA << '\n';
//    cout << "XA*XAi     =  \n"  << XA.contract(XAi, idxlist1{idx2(1, 0)}) << '\n';
//    cout << "YA         =  \n"  << YA << '\n';
//    cout << "YAi*YA     =  \n"  << YAi.contract(YA, idxlist1{idx2(1, 0)}) << '\n';
//
//    cout << "XB         =  \n"  << XB << '\n';
//    cout << "XB*XBi     =  \n"  << XB.contract(XBi, idxlist1{idx2(1, 0)}) << '\n';
//    cout << "YB         =  \n"  << YB << '\n';
//    cout << "YBi*YB     =  \n"  << YBi.contract(YB, idxlist1{idx2(1, 0)}) << '\n';
//
//    cout << "XA dims: " << XA.dimensions() << endl;
//    cout << "YA dims: " << YA.dimensions() << endl;
//    cout << "XB dims: " << XB.dimensions() << endl;
//    cout << "YB dims: " << YB.dimensions() << endl;

    //Symmetric
    Tensor2d YA_LA_XB = YA
            .contract(asDiagonal(MPS.LA), idxlist1{idx2(1, 0)})
            .contract(XB, idxlist1{idx2(1, 0)});
    //Symmetric
    Tensor2d YB_LB_XA = YB
            .contract(asDiagonal(MPS.LB), idxlist1{idx2(1, 0)})
            .contract(XA, idxlist1{idx2(1, 0)});

    SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
    chiA        = std::min(SVD.rank(),chi_max);
    Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
    Tensor2d VA = MatrixXd_to_Tensor2d(SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, U and V are not symmetric!
    Tensor1d LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());

    SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
    chiB        = std::min(SVD.rank(),chi_max);
    Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
    Tensor2d VB = MatrixXd_to_Tensor2d(SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, U and V are not symmetric!
    Tensor1d LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());
    cout << "LA new: \n" << LA << endl;
    cout << "LB new: \n" << LB << endl;
    cout << "S_E A: " <<  -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
    cout << "S_E B: " <<  -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;

//    cout << "UA: \n "<< UA << endl;
//    cout << "VA: \n "<< VA << endl;
//    cout << "UB: \n "<< UB << endl;
//    cout << "VB: \n "<< VB << endl;
    Tensor3d GA = VB
               .contract(XAi, idxlist1{idx2(1,0)})
               .contract(MPS.GA, idxlist1{idx2(1,1)})
               .contract(YAi, idxlist1{idx2(2, 0)})
               .contract(UA, idxlist1{idx2(2,0)})
               .shuffle(array3{1,0,2});
    Tensor3d GB = VA
            .contract(XBi, idxlist1{idx2(1,0)})
            .contract(MPS.GB, idxlist1{idx2(1,1)})
            .contract(YBi, idxlist1{idx2(2, 0)})
            .contract(UB, idxlist1{idx2(2,0)})
            .shuffle(array3{1,0,2});

    GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
    LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
    GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
    LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
    TA_R = GALA.contract(GALA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    TA_L = LBGA.contract(LBGA.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    TB_R = GBLB.contract(GBLB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});
    TB_L = LAGB.contract(LAGB.conjugate(), idxlist1{idx2(0, 0)}).shuffle(array4{0,2,1,3});


    VA_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TA_R.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiB, chiB});
    VA_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TA_L.reshape(array2{chiB*chiB,chiA*chiA})).real().reshape(array2{chiA, chiA});
    VB_R = Math::dominant_eigenvector<Math::Handedness::RGHT>  (TB_R.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiA, chiA});
    VB_L = Math::dominant_eigenvector<Math::Handedness::LEFT>  (TB_L.reshape(array2{chiA*chiA,chiB*chiB})).real().reshape(array2{chiB, chiB});


    cout << "VA_R: \n " << VA_R << '\n';
    cout << "VA_L: \n " << VA_L << '\n';
    cout << "VB_R: \n " << VB_R << '\n';
    cout << "VB_L: \n " << VB_L << '\n';
    //Check full 2-site Transfer matrix
    MatrixXd IL = MatrixXd::Identity(TA_L.dimension(0), TA_L.dimension(1));
    MatrixXd IR = MatrixXd::Identity(TB_R.dimension(2), TB_R.dimension(3));
    Tensor2d TL = TA_L.contract(TB_L, idxlist2{idx2(2,0), idx2(3,1)}).contract(MatrixXd_to_Tensor2d(IL), idxlist2{idx2(0,0),idx2(1,1)});
    Tensor2d TR = TA_R.contract(TB_R, idxlist2{idx2(2,0), idx2(3,1)}).contract(MatrixXd_to_Tensor2d(IR), idxlist2{idx2(2,0),idx2(3,1)});

    cout << "TL: \n" << TL << endl;
    cout << "TR: \n" << TR << endl;
//    MPS.GA = GA;
//    MPS.GB = GB;
//    MPS.LA = LA;
//    MPS.LB = LB;
//    MPS.L_tail = LB;

}



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
     cout << setprecision(16) << endl;
     cout << "LA old: \n" << MPS.LA << endl;
     cout << "LB old: \n" << MPS.LB << endl;
     cout << "S_E A: " <<  -MPS.LA.square().contract(MPS.LA.square().log().eval(), idxlist1{idx2(0,0)}) << endl;
     cout << "S_E B: " <<  -MPS.LB.square().contract(MPS.LB.square().log().eval(), idxlist1{idx2(0,0)}) << endl;

     Eigen::BDCSVD<MatrixXd> SVD;
//    SelfAdjointEigenSolver<MatrixXd> es;
//
     SVD.setThreshold(1e-24);
     int iter = 0;
     Tensor3d GA = MPS.GA;
     Tensor3d GB = MPS.GB;
     Tensor1d LA = MPS.LA;
     Tensor1d LB = MPS.LB;
     update_bond_dimensions();

     while (diffA + diffB > 1e-10) {
//        cout << "diff: " << diff << endl;
         iter++;
         //Coarse grain
         Tensor3d GA_old = GA;
         Tensor3d GB_old = GB;
         Tensor1d LA_old = LA;
         Tensor1d LB_old = LB;

         long chiA = LA.size();
         long chiB = LB.size();
         Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
         Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
         Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
         Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)

//        cout << "GALA dims: " << GALA.dimensions() << endl;
//        cout << "LBGA dims: " << LBGA.dimensions() << endl;
//        cout << "GBLB dims: " << GBLB.dimensions() << endl;
//        cout << "LAGB dims: " << LAGB.dimensions() << endl;


         Tensor2d VA_R = GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)});
         Tensor2d VA_L = LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)});
         Tensor2d VB_R = GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)});
         Tensor2d VB_L = LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)});

//        cout << "VA_R dims: " << VA_R.dimensions() << endl;
//        cout << "VA_L dims: " << VA_L.dimensions() << endl;
//        cout << "VB_R dims: " << VB_R.dimensions() << endl;
//        cout << "VB_L dims: " << VB_L.dimensions() << endl;
//         cout << "VA_R: \n " << VA_R << '\n';
//         cout << "VA_L: \n " << VA_L << '\n';
//         cout << "VB_R: \n " << VB_R << '\n';
//         cout << "VB_L: \n " << VB_L << '\n';



         //============================================================================//
         // Find square root using SVD decomposition
         //============================================================================//
//         SVD.compute(Tensor2d_to_MatrixXd(VA_R), Eigen::ComputeThinU);
//         Tensor2d XA = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//         Tensor2d XAi = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() *
//                 SVD.matrixU().transpose());
//         SVD.compute(Tensor2d_to_MatrixXd(VA_L), Eigen::ComputeThinU);
//         Tensor2d YA = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//         Tensor2d YAi = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() *
//                 SVD.matrixU().transpose());
//         SVD.compute(Tensor2d_to_MatrixXd(VB_R), Eigen::ComputeThinU);
//         Tensor2d XB = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//         Tensor2d XBi = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() *
//                 SVD.matrixU().transpose());
//         SVD.compute(Tensor2d_to_MatrixXd(VB_L), Eigen::ComputeThinU);
//         Tensor2d YB = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseSqrt().asDiagonal() * SVD.matrixU().transpose());
//         Tensor2d YBi = MatrixXd_to_Tensor2d(
//                 SVD.matrixU() * SVD.singularValues().cwiseInverse().cwiseSqrt().asDiagonal() *
//                 SVD.matrixU().transpose());

         //============================================================================//
         // Find square root using Eigenvalue decomposition
         //============================================================================//
         SelfAdjointEigenSolver<MatrixXd> es;
         es.compute(Tensor2d_to_MatrixXd(VA_R));
         Tensor2d XA  = MatrixXd_to_Tensor2d(es.operatorSqrt());
         Tensor2d XAi = MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
         es.compute(Tensor2d_to_MatrixXd(VA_L));
         Tensor2d YA  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
         Tensor2d YAi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
         es.compute(Tensor2d_to_MatrixXd(VB_R));
         Tensor2d XB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
         Tensor2d XBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());
         es.compute(Tensor2d_to_MatrixXd(VB_L));
         Tensor2d YB  =  MatrixXd_to_Tensor2d(es.operatorSqrt());
         Tensor2d YBi =  MatrixXd_to_Tensor2d(es.operatorInverseSqrt());

         //Symmetric
         Tensor2d YA_LA_XB = YA
                 .contract(asDiagonal(LA), idxlist1{idx2(1, 0)})
                 .contract(XB, idxlist1{idx2(1, 0)});
         //Symmetric
         Tensor2d YB_LB_XA = YB
                 .contract(asDiagonal(LB), idxlist1{idx2(1, 0)})
                 .contract(XA, idxlist1{idx2(1, 0)});

         SVD.compute(Tensor2d_to_MatrixXd(YA_LA_XB), Eigen::ComputeThinU | Eigen::ComputeThinV);
         chiA = std::min(SVD.rank(), chi_max);
         Tensor2d UA = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiA));
         Tensor2d VA = MatrixXd_to_Tensor2d(
                 SVD.matrixV().leftCols(chiA).transpose()); //Transpose is important, U and V are not symmetric!
         LA = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiA).normalized());

         SVD.compute(Tensor2d_to_MatrixXd(YB_LB_XA), Eigen::ComputeThinU | Eigen::ComputeThinV);
         chiB = std::min(SVD.rank(), chi_max);
         Tensor2d UB = MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chiB));
         Tensor2d VB = MatrixXd_to_Tensor2d(
                 SVD.matrixV().leftCols(chiB).transpose()); //Transpose is important, U and V are not symmetric!
         LB = MatrixXd_to_Tensor1d(SVD.singularValues().head(chiB).normalized());

//    cout << "UA: \n "<< UA << endl;
//    cout << "VA: \n "<< VA << endl;
//    cout << "UB: \n "<< UB << endl;
//    cout << "VB: \n "<< VB << endl;
         GA = VB
                 .contract(XAi, idxlist1{idx2(1, 0)})
                 .contract(GA_old, idxlist1{idx2(1, 1)})
                 .contract(YAi, idxlist1{idx2(2, 0)})
                 .contract(UA, idxlist1{idx2(2, 0)})
                 .shuffle(array3{1, 0, 2});
         GB = VA
                 .contract(XBi, idxlist1{idx2(1, 0)})
                 .contract(GB_old, idxlist1{idx2(1, 1)})
                 .contract(YBi, idxlist1{idx2(2, 0)})
                 .contract(UB, idxlist1{idx2(2, 0)})
                 .shuffle(array3{1, 0, 2});


         Tensor0d tempA = (LA - LA_old.slice(array1{0}, array1{chiA})).square().sum().sqrt();
         Tensor0d tempB = (LB - LB_old.slice(array1{0}, array1{chiB})).square().sum().sqrt();
         diffA = tempA(0);
         diffB = tempB(0);
//         cout << "LA diff: " << diffA << endl;
//         cout << "LB diff: " << diffB << endl;
//        if(iter > 100 ){
//            exit(0);
//        }
     }
     Tensor3d GALA = GA.contract(asDiagonal(LA), idxlist1{idx2(2, 0)});                          //(eq 107)
     Tensor3d LBGA = asDiagonal(LB).contract(GA, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 106)
     Tensor3d GBLB = GB.contract(asDiagonal(LB), idxlist1{idx2(2, 0)});                          //(eq 109)
     Tensor3d LAGB = asDiagonal(LA).contract(GB, idxlist1{idx2(1, 1)}).shuffle(array3{1, 0, 2}); //(eq 108)
     cout << "TA_R: \n " << GALA.contract(GALA.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
     cout << "TA_L: \n " << LBGA.contract(LBGA.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
     cout << "TB_R: \n " << GBLB.contract(GBLB.conjugate(), idxlist2{idx2(0, 0), idx2(2, 2)}) << '\n';
     cout << "TB_L: \n " << LAGB.contract(LAGB.conjugate(), idxlist2{idx2(0, 0), idx2(1, 1)}) << '\n';
//    Tensor2d LGL =   asDiagonal(L)
//                   .contract(G, idxlist1{idx2(1,2)})
//                   .contract(asDiagonal(L), idxlist1{idx2(3,0)})
//                   .shuffle(array4{1,0,2,3})
//                   .reshape(array2{G.dimension(0)*G.dimension(2), G.dimension(1)*G.dimension(3)});
//
//
////    //Note the difference! U3 is a d x [chib] x chi matrix! Previously we had [chia] in normal SVD truncation
//    SVD.compute(tensor2_to_matrix(LGL), Eigen::ComputeThinU | Eigen::ComputeThinV);
//    chiL                         = std::min(SVD.rank(),chi_max);
//    truncation_error             = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiL).sum();
//    Tensor3d U3                   = matrix_to_tensor<3>(SVD.matrixU(),{d,chiR,chiL});
//    MPS.LA                       = matrix_to_tensor<1>(SVD.singularValues().head(chiL).normalized(), {chiL});
//    Tensor3d V3                   = matrix_to_tensor<3>(SVD.matrixV().leftCols(chiL),{d,chiR,chiL}).shuffle(array3{0,2,1});
////
//    MPS.GA = asInverseDiagonal(L).contract(U3, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
//    MPS.GB = V3.contract(asInverseDiagonal(L), idxlist1{idx2(2,0)});
//    MPS.LB = L;
//    MPS.L_tail = L;

//    cout << "Full LA after: \n" << SVD.singularValues().normalized() << endl;
     cout << "LA after: \n" << LA << endl;
     cout << "LB after: \n" << LB << endl;
     cout << "Entropy A: " << -LA.square().contract(LA.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
     cout << "Entropy B: " << -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0, 0)}) << endl;
     cout << "\n\niter: " << iter << endl;

//     MPS.GA = GA;
//     MPS.GB = GB;
//     MPS.LA = LA;
//     MPS.LB = LB;
//     MPS.L_tail = LB;
//    TRY USING ROBUST SVD!
//             Stoudenmire, E. M., & White, S. R. (2013). Real-space parallel density matrix renormalization group. Physical Review B, 87(April), 155137. https://doi.org/10.1103/PhysRevB.87.155137
}
