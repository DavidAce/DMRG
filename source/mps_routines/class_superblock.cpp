//
// Created by david on 7/22/17.
//
//#define EIGEN_USE_MKL_ALL

#include "class_superblock.h"
#include "class_eigensolver_product.h"
#include <iomanip>
#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <sim_parameters/n_model.h>
#include <general/class_svd_wrapper.h>



using namespace Textra;
using namespace Eigen;


class_superblock::class_superblock(){
    MPS.initialize(H.local_dimension);
    update_bond_dimensions();
}

//============================================================================//
// Find smallest eigenvalue using Spectra.
//============================================================================//
void class_superblock::find_ground_state(int eigSteps, double eigThreshold){
    Textra::Tensor<1,double> theta = MPS.get_theta().real().shuffle(array4{0,2,1,3}).reshape(shape1);
    class_eigensolver_product esp(Lblock.block, Rblock.block, H.MM, shape4);
    Spectra::SymEigsSolver<double,
            Spectra::SMALLEST_ALGE,
            class_eigensolver_product>
            eigs(&esp, 1, std::min(10,(int)theta.size()) );

//    Spectra::SymEigsSolver<double,
//            Spectra::SMALLEST_ALGE,
//            class_superblock>
//            eigs(this, 1, std::min(10,(int)theta.size()) );
//    cout << setprecision(16) << "theta:\n" << MPS.get_theta().shuffle(array4{0,2,1,3}).reshape(shape2) << endl;
    eigs.init(theta.data()); //Throw in the initial vector v0 here
    eigs.compute(eigSteps,eigThreshold, Spectra::SMALLEST_ALGE);

    if(eigs.info() != Spectra::SUCCESSFUL){
        cout << "Eigenvalue solver failed." << '\n';
        exit(1);
    }
    ground_state =  Matrix_to_Tensor<2,Scalar>(eigs.eigenvectors(), shape2);

}

//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS.
//============================================================================//
void class_superblock::truncate(long chi_max_, double SVDThreshold){
    class_SVD<Scalar> SVD;
    SVD.setThreshold(SVDThreshold);
    chi_max = chi_max_;
    std::tie(U, MPS.LA, V) = SVD.schmidt(ground_state, d, chiL, chi_max, chiR);
    truncation_error = SVD.get_truncation_error();
}




void class_superblock::time_evolve(){
//    ground_state  = H.asTimeEvolution.contract(MPS.get_theta(), idxlist2 {idx2(0,1),idx2(1,2)}).shuffle(array4{0,2,1,3}).reshape(shape2);
    ground_state  = H.Udt.contract(MPS.get_theta(), idx<2>({0,1},{0,1}))
                                   .shuffle(array4{0,2,1,3})
                                   .reshape(shape2);
}


void class_superblock::update_MPS(){
    MPS.GA  = asInverseDiagonal(MPS.L_tail).contract(U,idx<1>({1},{1})).shuffle(array3{1,0,2});
    MPS.GB  = V.contract(asInverseDiagonal(MPS.LB), idx<1>({2},{0}));
}

void class_superblock::enlarge_environment(int direction){
    if (direction == 1){
        Lblock.enlarge(MPS, H.M);
    }else if (direction == -1){
        Rblock.enlarge(MPS, H.M);
    }else if(direction == 0){
        Lblock.enlarge(MPS, H.M);
        Rblock.enlarge(MPS, H.M);
        chain_length += 2;
    }
}


void class_superblock::initialize() {
    MPS.GA.resize(array3{d,1,1});
    MPS.GA.setZero();
    MPS.GA(0,0,0) = 0;
    MPS.GA(1,0,0) = 1;
    MPS.LA.resize(array1{1});
    MPS.LA.setConstant(1.0);

    MPS.GB.resize(array3{d,1,1});
    MPS.GB.setZero();
    MPS.GB(0,0,0) = 0;
    MPS.GB(1,0,0) = 1;
    MPS.LB.resize(array1{1});
    MPS.LB.setConstant(1.0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
    enlarge_environment(0);
}


void class_superblock::update_bond_dimensions(){
    d       = H.local_dimension;
    chi     = MPS.GA.dimension(2);
    chiL    = MPS.GA.dimension(1);
    chiR    = MPS.GB.dimension(2);
    shape1  = {d * chiL * d * chiR};
    shape2  = {d * chiL , d * chiR};
    shape4  = {d , chiL , d , chiR};
}


void  class_superblock::swap_AB(){
    MPS.swap_AB();
}

void class_superblock::reset(){
    *this = class_superblock();
}



//Functions for eigenvalue solver Spectra
//int class_superblock::rows()const {return (int)shape1[0];}
//int class_superblock::cols()const {return (int)shape1[0];}
//void class_superblock::perform_op(const double *x_in, double *y_out) const {
//    Eigen::TensorMap<Textra::const_Tensor<4,double>> x_input (x_in, shape4);
//    Textra::Tensor<1,double> y_output = Lblock.block.contract(x_input.cast<Scalar>(),     idx<1>({1},{1}))
//                                    .contract(H.WW ,       idx<3>({1,2,3},{0,4,5}))
//                                    .contract(Rblock.block,idx<2>({1,2},{1,2}))
//                                    .shuffle(array4{1,0,2,3})
//                                    .reshape(shape1).real();
//    std::move(y_output.data(), y_output.data()+y_output.size(),  y_out);
//}
