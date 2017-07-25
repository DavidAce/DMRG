//
// Created by david on 7/22/17.
//

#include "class_superblock.h"
class_superblock::class_superblock(const int eigSteps_,
                 const double eigThreshold_,
                 const double SVDThreshold_,
                 const long chi_){
    MPS.initialize(H.local_dimension, H.sites);
    eigSteps = eigSteps_;
    eigThreshold = eigThreshold_;
    SVDThreshold = SVDThreshold_;
    chi = chi_;
}

void class_superblock::find_ground_state(){
    //============================================================================//
    // Find smallest eigenvalue using Spectra
    //============================================================================//
    Spectra::SymEigsSolver<double,
            Spectra::SMALLEST_ALGE,
            class_superblock>
            eigs(this, 1, 4);
    array1 tempshape = {MPS.G.A.dimension(1) * MPS.G.B.dimension(2) * H.local_dimension * H.local_dimension};
    cout << "Shape1 theta: " << tempshape[0] << endl;

    Tensor1 theta = MPS.get_theta().reshape(shape1);
    eigs.init(theta.data()); //Throw in the initial vector v0 here
    eigs.compute(eigSteps,eigThreshold, Spectra::SMALLEST_ALGE);

    if(eigs.info() != Spectra::SUCCESSFUL){
        cout << "Eigenvalue solver failed." << endl;
        exit(1);
    }
    ground_state =  matrix_to_tensor2(eigs.eigenvectors(), shape2);
}

void class_superblock::truncate(){
    std::tie(U, MPS.L.A, V) = SVD.decompose(ground_state, SVDThreshold, H.local_dimension, chi, chia, chib);
}


void class_superblock::update_MPS(){
    MPS.G.A  = asInverseDiagonal(MPS.L_tail).contract(U, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
    MPS.G.B  = V.contract(asInverseDiagonal(MPS.L.B), idxlist1{idx2(2,0)});
}

void class_superblock::enlarge_environment(){
    Lblock.enlarge(MPS, H.W);
    Rblock.enlarge(MPS, H.W);
}

void class_superblock::enlarge_environment(long direction){
    if (direction == 1){
        Lblock.enlarge(MPS, H.W);
    }else{
        Rblock.enlarge(MPS, H.W);
    }
}

void class_superblock::update_bond_dimensions(){
    chia = MPS.G.A.dimension(1);
    chib = MPS.G.B.dimension(2);
    shape1 = {H.local_dimension * chia * H.local_dimension * chib};
    shape2 = {H.local_dimension * chia , H.local_dimension * chib};
    shape4 = {H.local_dimension , chia , H.local_dimension , chib};
}


void  class_superblock::swap_AB(){
    MPS.swap_AB();
}

void class_superblock::print_picture(){
    cout << Lblock.full_picture << H.pic << Rblock.full_picture << endl;
}

void class_superblock::print_error() {
    MPS.print_error(H.asTensor4);
}


//Functions for eigenvalue solver Spectra
int class_superblock::rows()const {return (int)shape1[0];}
int class_superblock::cols()const {return (int)shape1[0];}

void class_superblock::perform_op(const double *x_in, double *y_out) const {
    Eigen::TensorMap<const_Tensor4> x_input (x_in, shape4);
    Tensor1 y_output = Lblock.block.contract(x_input, idxlist1{idx2(0,1)})
            .contract(H.W , idxlist2{idx2(1,0),idx2(2,2)})
            .contract(H.W , idxlist2{idx2(3,0),idx2(1,2)})
            .contract(Rblock.block, idxlist2{idx2(1,0),idx2(3,2)})
            .shuffle(array4{1,0,2,3})
            .reshape(shape1);
    std::move(y_output.data(), y_output.data()+y_output.size(),  y_out);
}
