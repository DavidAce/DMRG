//
// Created by david on 7/22/17.
//

#include <class_superblock.h>
#include <iomanip>
#include "n_model.h"


using namespace Textra;

class_superblock::class_superblock( const int eigSteps_,
                                    const double eigThreshold_,
                                    const double SVDThreshold_,
                                    const long chi_){
    MPS.initialize(H.local_dimension);
    chi             = chi_;
    update_bond_dimensions();
    eigSteps        = eigSteps_;
    eigThreshold    = eigThreshold_;
    SVDThreshold    = SVDThreshold_;

}

void class_superblock::find_ground_state(){
    //============================================================================//
    // Find smallest eigenvalue using Spectra
    //============================================================================//
    Spectra::SymEigsSolver<double,
            Spectra::SMALLEST_ALGE,
            class_superblock>
            eigs(this, 1, 4);

    Textra::Tensor1 theta = MPS.get_theta().reshape(shape1);
    eigs.init(theta.data()); //Throw in the initial vector v0 here
    eigs.compute(eigSteps,eigThreshold, Spectra::SMALLEST_ALGE);

    if(eigs.info() != Spectra::SUCCESSFUL){
        cout << "Eigenvalue solver failed." << endl;
        exit(1);
    }
    ground_state =  matrix_to_tensor2(eigs.eigenvectors(), shape2);
}

void class_superblock::time_evolve(){
    ground_state  = H.asTimeEvolution.contract(MPS.get_theta(), idxlist2 {idx2(2,0),idx2(3,2)}).shuffle(array4{0,2,1,3}).reshape(shape2);

}

void class_superblock::truncate(){
    Eigen::BDCSVD<MatrixType> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(tensor2_to_matrix(ground_state), Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi2       = std::min(SVD.rank(),chi);

//    cout << "chi chi2 chia chib: " << chi  << " " << chi2 << " " << chia << " " << chib << endl;
    double renorm   = SVD.singularValues().head(chi2).norm();
    U               = matrix_to_tensor3(SVD.matrixU().leftCols(chi2),{d,chia,chi2});
    MPS.LA          = matrix_to_tensor1(SVD.singularValues().head(chi2)) / renorm;
    V               = matrix_to_tensor3(SVD.matrixV().leftCols(chi2),{d,chib,chi2}).shuffle(array3{0,2,1});
}


void class_superblock::update_MPS(){
    MPS.GA  = asInverseDiagonal(MPS.L_tail).contract(U, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
    MPS.GB  = V.contract(asInverseDiagonal(MPS.LB), idxlist1{idx2(2,0)});
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
    d       = H.local_dimension;
    chia    = MPS.GA.dimension(1);
    chib    = MPS.GB.dimension(2);
    shape1  = {H.local_dimension * chia * H.local_dimension * chib};
    shape2  = {H.local_dimension * chia , H.local_dimension * chib};
    shape4  = {H.local_dimension , chia , H.local_dimension , chib};
}


void  class_superblock::swap_AB(){
    MPS.swap_AB();
}

void class_superblock::reset(){
    *this = class_superblock(eigSteps, eigThreshold, SVDThreshold, chi);
}

void class_superblock::print_picture(){
    cout << Lblock.picture << H.picture << Rblock.picture << endl;
}

void class_superblock::print_error_DMRG() {
    double energy = MPS.get_energy(H.asTensor4);
    cout << setprecision(16);
    cout << "E_iDMRG = " << energy << endl;
    cout << "E_exact = " <<  Model::energy_exact << endl;
    cout << std::scientific;
    cout << "Error   = " <<  Model::energy_exact - energy << endl;
}



void class_superblock::print_error_TEBD() {
    double energy = MPS.get_energy(H.asTensor4);
    cout << setprecision(16);
    cout << "E_iTEBD = " << energy << endl;
    cout << "E_exact = " << Model::energy_exact << endl;
    cout << std::scientific;
    cout << "Error   = " << Model::energy_exact - energy << endl;
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
