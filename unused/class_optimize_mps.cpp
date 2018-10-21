//
// Created by david on 2018-05-04.
//

#include "class_optimize_mps.h"
#include <arpackpp/arscomp.h>
#include <sim_parameters/nmspc_sim_settings.h>
using Scalar = std::complex<double>;

void class_optimize_mps::optimize_mps(int nev)
{
    class_contraction<Scalar>  contraction(Lblock, Rblock, HA, HB, shape_theta4,shape_mpo4,testmat);
    int dim = contraction.cols();
    int size = contraction.cols() * contraction.rows();
    int ncv = std::min(settings::precision::eig_max_ncv,size);
//    std::cout << "dim: " << dim << " size: " << size << " ncv: " << ncv << " nev: " << nev << std::endl;
//    std::cout << "shape_theta4: " << shape_theta4[0]  << " " << shape_theta4[1]  << " " << shape_theta4[2]  << " " << shape_theta4[3] << std::endl;
//    std::cout << "shape_mpo4  : " << shape_mpo4[0]  << " " << shape_mpo4[1]  << " " << shape_mpo4[2]  << " " << shape_mpo4[3] << std::endl;
//    std::cout << "shape_mpo4" << shape_mpo4 << std:endl;

    ARCompStdEig<double, class_contraction<Scalar>> eig (dim, nev, &contraction, &class_contraction<Scalar>::MultMv, "SR", ncv,settings::precision::eigThreshold,settings::precision::eigMaxIter, resid);
    eig.FindEigenvectors();
    rows = eig.GetN();
    cols = eig.GetNev();
    iter = eig.GetIter();
    counter = contraction.counter;
    eigvals = std::vector<Scalar>(eig.RawEigenvalues() , eig.RawEigenvalues()  + cols);
    eigvecs = std::vector<Scalar>(eig.RawEigenvector(0), eig.RawEigenvectors() + rows*cols);
}