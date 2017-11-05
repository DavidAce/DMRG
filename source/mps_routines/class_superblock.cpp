//
// Created by david on 7/22/17.
//
#include "class_superblock.h"
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
    class_SVD SVD;
    chi_max = chi_max_;
    std::tie(U, MPS.LA, V) = SVD.schmidt(ground_state, d, chiL, chi_max, chiR);
    truncation_error = SVD.get_truncation_error();
}




void class_superblock::time_evolve(){
//    ground_state  = H.asTimeEvolution.contract(MPS.get_theta(), idxlist2 {idx2(0,1),idx2(1,2)}).shuffle(array4{0,2,1,3}).reshape(shape2);
    ground_state  = MPS.get_theta().contract( H.asTimeEvolution, idxlist2 {idx2(1,0),idx2(2,1)}).shuffle(array4{2,0,3,1}).reshape(shape2);
}


void class_superblock::update_MPS(){
    MPS.GA  = asInverseDiagonal(MPS.L_tail).contract(U, idxlist1{idx2(1,1)}).shuffle(array3{1,0,2});
    MPS.GB  = V.contract(asInverseDiagonal(MPS.LB), idxlist1{idx2(2,0)});
}

void class_superblock::enlarge_environment(int direction){
    if (direction == 1){
        Lblock.enlarge(MPS, H.W);
    }else if (direction == -1){
        Rblock.enlarge(MPS, H.W);
    }else if(direction == 0){
        Lblock.enlarge(MPS, H.W);
        Rblock.enlarge(MPS, H.W);
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



double class_superblock::get_expectationvalue(const Tensor4d &MPO){
    Tensor4d theta = MPS.get_theta();
    Tensor0d result = Lblock.block
            .contract(theta.conjugate(), idxlist1{idx2(0,0)})
            .contract(MPO,   idxlist2{idx2(1,0), idx2(2,2)})
            .contract(MPO,   idxlist2{idx2(3,0), idx2(1,2)})
            .contract(theta, idxlist3{idx2(0,0), idx2(2,1), idx2(4,2)})
            .contract(Rblock.block, idxlist3{idx2(0,0), idx2(1,2), idx2(2,1)});
    return result(0)/chain_length;
}

double class_superblock::get_energy(){
    return MPS.get_energy(H.asTensor4);
}

double class_superblock::get_entropy(){
    return MPS.get_entropy();
}

double class_superblock::get_variance(){
//    return MPS.get_variance(H.asTensor4);
    return MPS.get_variance(H.asTensor4);
}


void class_superblock::update_bond_dimensions(){
    assert(MPS.GA.dimension(2) == MPS.LA.size() && "GA LA Dimension mismatch");
    assert(MPS.GB.dimension(1) == MPS.LA.size() && "GB LA Dimension mismatch");
    assert(MPS.GA.dimension(1) == MPS.LB.size() && "GA LB Dimension mismatch");
    assert(MPS.GB.dimension(2) == MPS.LB.size() && "GB LB Dimension mismatch");
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

void class_superblock::print_picture(bool graphics){
    if (graphics) {
        cout << Lblock.picture << H.picture << Rblock.picture << '\n';
    }
}

void class_superblock::print_status_full(int verbosity, string name) {
    if (verbosity >= 1) {
//        double energy_fin  = get_expectationvalue(H.W);
        double energy  = get_energy();
        double entropy  = get_entropy();
        double variance = get_variance();
        cout << setprecision(16) << fixed << left;
//        cout << setw(20) << "E_fin" + name << " = "   << energy_fin << '\n';
        cout << setw(20) << "E_inf_" + name << " = "   << energy << '\n';
        if(verbosity >= 2) {
            cout  << setw(20) << "E_exact"              << " = " << setprecision(16) << fixed      << Model::energy_exact << '\n';
            cout  << setw(20) << "Entanglement Entropy" << " = " << setprecision(16) << fixed      << entropy << '\n';
            cout  << setw(20) << "E_error"              << " = " << setprecision(4)  << scientific << energy - Model::energy_exact << '\n';
            cout  << setw(20) << "Truncation error"     << " = " << setprecision(4)  << scientific << truncation_error << '\n';
            cout  << setw(20) << "Variance"             << " = " << setprecision(16) << scientific << variance << '\n';
            cout  << setw(20) << "chi_max"              << " = " << setprecision(4)  << fixed      << chi_max << '\n';
            cout  << setw(20) << "chi    "              << " = " << setprecision(4)  << fixed      << chi     << '\n';

        }
    }
}


void class_superblock::print_status_update(string name, int step) {
    double energy  = get_energy();
    double entropy  = get_entropy();
    double variance = get_variance();
    cout << left  << name
         << left  << setw(10) << step
         << left  << "E: "              << setw(25) << setprecision(16)    << fixed   << energy
         << left  << "E-E_exact: "      << setw(25) << setprecision(16)    << fixed   << energy - Model::energy_exact
         << left  << "Trunc. error: "   << setw(25) << setprecision(16)    << fixed   << truncation_error
         << left  << "S_E: "            << setw(25) << setprecision(16)    << fixed   << entropy
         << left  << "Var: "            << setw(25) << setprecision(16)    << fixed   << variance
         << left  << "chi_max: "        << setw(3)  << setprecision(4)     << fixed   << chi_max
         << left  << "chi: "            << setw(20) << setprecision(4)     << fixed   << chi     << '\n';
}



//Functions for eigenvalue solver Spectra
int class_superblock::rows()const {return (int)shape1[0];}
int class_superblock::cols()const {return (int)shape1[0];}
void class_superblock::perform_op(const double *x_in, double *y_out) const {
    Eigen::TensorMap<const_Tensor4d> x_input (x_in, shape4);
    Tensor1d y_output = Lblock.block.contract(x_input, idxlist1{idx2(0,1)})
            .contract(H.W , idxlist2{idx2(1,0),idx2(2,2)})
            .contract(H.W , idxlist2{idx2(3,0),idx2(1,2)})
            .contract(Rblock.block, idxlist2{idx2(1,0),idx2(3,2)})
            .shuffle(array4{1,0,2,3})
            .reshape(shape1);

//    Tensor1d y_output = Lblock.block.contract(x_input, idxlist1{idx2(1,1)})
//            .contract(H.W , idxlist2{idx2(1,0),idx2(2,3)})
//            .contract(H.W , idxlist2{idx2(3,0),idx2(1,3)})
//            .contract(Rblock.block, idxlist2{idx2(1,1),idx2(3,2)})
//            .shuffle(array4{1,0,2,3})
//            .reshape(shape1);
    std::move(y_output.data(), y_output.data()+y_output.size(),  y_out);
}
