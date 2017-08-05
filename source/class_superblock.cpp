//
// Created by david on 7/22/17.
//
#include <class_superblock.h>
#include <iomanip>
#include <Eigen/SVD>
#include <Eigen/Core>
#include <n_model.h>


using namespace Textra;

class_superblock::class_superblock(){
    MPS.initialize(H.local_dimension);
    update_bond_dimensions();
}

void class_superblock::find_ground_state(int eigSteps, double eigThreshold){
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

void class_superblock::truncate(long chi, double SVDThreshold){
    Eigen::BDCSVD<MatrixType> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(tensor2_to_matrix(ground_state), Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi2       = std::min(SVD.rank(),chi);

//    cout << "chi chi2 nonzero rank singv: " << chi  << " " << chi2 << " " << SVD.nonzeroSingularValues() << " " << SVD.rank()<< " "<< SVD.singularValues().size() << endl;
//    double renorm   = 1.0 - SVD.singularValues().tail(SVD.singularValues().size()-chi2).sum();
//    truncation_error= SVD.singularValues().tail(SVD.singularValues().size()-chi2).sum();
    double renorm   = SVD.singularValues().head(chi2).norm();
    truncation_error= 1.0-renorm;
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
    *this = class_superblock();
}

void class_superblock::print_picture(bool graphics){
    if (graphics) {
        cout << Lblock.picture << H.picture << Rblock.picture << endl;
    }
}

void class_superblock::print_state(int verbosity, string name) {
    if (verbosity >= 1) {
        double energy   = MPS.get_energy(H.asTensor4);
        double entropy  = MPS.get_entropy();
        cout << setprecision(16) << fixed << left;
        cout << setw(20) << "E_" + name << " = "   << energy << endl;
        if(verbosity >= 2) {
            cout  << setw(20) << "E_exact"              << " = " << setprecision(16) << fixed      << Model::energy_exact << endl;
            cout  << setw(20) << "Entanglement Entropy" << " = " << setprecision(16) << fixed      << entropy << endl;
            cout  << setw(20) << "E_error"              << " = " << setprecision(4)  << scientific << energy - Model::energy_exact << endl;
            cout  << setw(20) << "Truncation error"     << " = " << setprecision(4)  << scientific << truncation_error << endl;


        }
    }
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
