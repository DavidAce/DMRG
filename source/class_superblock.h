//
// Created by david on 7/22/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_SUPERBLOCK_H
#define FINITE_DMRG_EIGEN_CLASS_SUPERBLOCK_H


#include <class_environment.h>
//#include <class_SVD.h>
#include <class_MPS.h>
#include <class_TwoSiteHamiltonian.h>

class class_superblock {
private:
    //Bond dimensions and tensor shapes needed by the eigensolver and SVD.
    long chia,chib;
    array1 shape1;
    array2 shape2;
    array4 shape4;

    //Variables that control eigensolver and SVD precision
    int     eigSteps;
    double  eigThreshold;
    long    chi;
    double  SVDThreshold;

    //Store temporaries for eigensolver and SVD.
    Tensor2 ground_state;   //Stores the ground state eigenvector
    Tensor3 U,V;            //For SVD decomposition
public:

    class_environment_L Lblock;
    class_environment_R Rblock;
    class_MPS MPS;
    class_TwoSiteHamiltonian H;

    class_superblock(const int eigSteps_,
                     const double eigThreshold_,
                     const double SVDThreshold_,
                     const long chi_);

    void find_ground_state();
    void truncate();
    void update_MPS();
    void enlarge_environment();
    void enlarge_environment(long direction);
    void update_bond_dimensions();
    void swap_AB();
    void print_picture();
    void print_error();

    //Functions for eigenvalue solver Spectra
    void perform_op(const double *x_in, double *y_out) const;
    int rows()const;
    int cols()const;



};


#endif //FINITE_DMRG_EIGEN_CLASS_SUPERBLOCK_H
