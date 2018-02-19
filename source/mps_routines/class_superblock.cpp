//
// Created by david on 7/22/17.
//

#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <iomanip>
#include <SymEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <general/class_svd_wrapper.h>
#include <mps_routines/class_custom_contraction.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include "arssym.h"
#include "symsol.h"

using namespace std;
using namespace Textra;


class_superblock::class_superblock():
        MPS(std::make_shared<class_mps<Scalar>>()),
        H(std::make_shared<class_mpo<Scalar>>()),
        Lblock(std::make_shared<class_environment<Scalar>>("L")),
        Rblock(std::make_shared<class_environment<Scalar>>("R")),
        Lblock2(std::make_shared<class_environment_var<Scalar>>("L")),
        Rblock2(std::make_shared<class_environment_var<Scalar>>("R"))
{
    MPS->initialize(H->local_dimension);
    chain_length = H->mps_sites;
    update_bond_dimensions();
    initialize();
}

////============================================================================//
//// Find smallest eigenvalue using Spectra.
////============================================================================//
//void class_superblock::find_ground_state(int eigSteps, double eigThreshold){
//    Tensor<double,1> theta = MPS->get_theta().shuffle(array4{0,2,1,3}).reshape(shape1);
//    class_eigensolver_product  esp(theta,Lblock->block, Rblock->block, H->MM, shape4);
//    Spectra::SymEigsSolver<double,
//                           Spectra::SMALLEST_ALGE,
//                           class_eigensolver_product>
//                           eigs(&esp, 1, std::min(10,(int)theta.size()) );
//    eigs.init(theta.data()); //Throw in the initial vector v0 here
//    eigs.compute(eigSteps,eigThreshold, Spectra::SMALLEST_ALGE);
//
//    if(eigs.info() != Spectra::SUCCESSFUL){
//        cout << "Eigenvalue solver failed." << '\n';
//        exit(1);
//    }
//    ground_state =  Matrix_to_Tensor<Scalar,2>(eigs.eigenvectors(), shape2);
//
//}


//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//
void class_superblock::find_ground_state_arpack(int eigSteps, double eigThreshold){
    Tensor<Scalar,1> theta = MPS->get_theta().shuffle(array4{0,2,1,3}).reshape(shape1);
    class_custom_contraction<Scalar>  esp(Lblock->block, Rblock->block, H->MM, shape4);
    int ncv = std::min(settings::precision::eig_max_ncv,(int)theta.size());
    int nev = 1;
    int dim = esp.cols();
    ARSymStdEig<Scalar, class_custom_contraction<Scalar>> eig (dim, nev, &esp, &class_custom_contraction<Scalar>::MultMv, "SA", ncv,eigThreshold,eigSteps,theta.data() );
    eig.FindEigenvectors();
    using T_vec = decltype(eig.Eigenvector(0, 0));
    using T_val = decltype(eig.Eigenvalue(0));
    int rows = eig.GetN();
    int cols = std::min(1, eig.GetNev());
    Textra::Tensor<T_vec, 2> eigvecs(rows, cols);
    Textra::Tensor<T_val, 1> eigvals(cols);
    for (int i = 0; i < cols; i++) {
        eigvals(i) = eig.Eigenvalue(0);
        for (int j = 0; j < rows; j++) {
            eigvecs(j, i) = eig.Eigenvector(i, j);
        }
    }
    energy = eigvals(0);
    ground_state = eigvecs.reshape(shape2);
}


//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS->
//============================================================================//
void class_superblock::truncate(long chi_max_, double SVDThreshold){
    class_SVD<Scalar> SVD;
    SVD.setThreshold(SVDThreshold);
    std::tie(U, MPS->LA, V) = SVD.schmidt(ground_state, d, chiL, chi_max_, chiR);
    truncation_error = SVD.get_truncation_error();
}



//============================================================================//
// Do iTEBD time evolution
//============================================================================//
void class_superblock::time_evolve(){
    ground_state  = H->Udt.contract(MPS->get_theta(), idx<2>({0,1},{0,1}))
                                   .shuffle(array4{0,2,1,3})
                                   .reshape(shape2);
}


void class_superblock::update_MPS(){
    MPS->GA  = asDiagonalInversed(MPS->L_tail).contract(U,idx<1>({1},{1})).shuffle(array3{1,0,2});
    MPS->GB  = V.contract(asDiagonalInversed(MPS->LB), idx<1>({2},{0}));
}

void class_superblock::enlarge_environment(int direction){
    if (direction == 1){
        Lblock->enlarge(MPS, H->M);
        Lblock2->enlarge(MPS, H->M);
    }else if (direction == -1){
        Rblock->enlarge(MPS, H->M);
        Rblock2->enlarge(MPS, H->M);
    }else if(direction == 0){
        Lblock->enlarge(MPS, H->M);
        Rblock->enlarge(MPS, H->M);
        Lblock2->enlarge(MPS, H->M);
        Rblock2->enlarge(MPS, H->M);
        chain_length += 2;
    }
}


void class_superblock::initialize() {
    MPS->GA.resize(array3{d,1,1});
    MPS->GA.setZero();
    MPS->GA(1,0,0) = 1;
    MPS->GA(0,0,0) = 0;
    MPS->LA.resize(array1{1});
    MPS->LA.setConstant(1.0);

    MPS->GB.resize(array3{d,1,1});
    MPS->GB.setZero();
    MPS->GB(1,0,0) = 1;
    MPS->GB(0,0,0) = 0;
    MPS->LB.resize(array1{1});
    MPS->LB.setConstant(1.0);
    MPS->L_tail.resize(array1{1});
    MPS->L_tail.setConstant(1.0);
}





void class_superblock::update_bond_dimensions(){
    d       = H->local_dimension;
    chi     = MPS->GA.dimension(2);
    chiL    = MPS->GA.dimension(1);
    chiR    = MPS->GB.dimension(2);
    shape1  = {d * chiL * d * chiR};
    shape2  = {d * chiL , d * chiR};
    shape4  = {d , chiL , d , chiR};
}


void  class_superblock::swap_AB(){
    MPS->swap_AB();
}

//void class_superblock::reset(){
//    *this = class_superblock();
//}
