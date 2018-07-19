//
// Created by david on 2018-05-06.
//

#include "class_arpack_eigsolver.h"
#include <assert.h>
#include <algorithm>
#include <arpack++/ardsnsym.h>
#include <arpack++/ardscomp.h>
#include <arpack++/ardgcomp.h>
#include <arpack++/ardssym.h>
#include <arpack++/arseig.h>
#include <general/class_arpack_custom_products.h>
//#include <sim_parameters/nmspc_sim_settings.h>

#define profile_eigsolver 0


template<typename Scalar, Form form>
class_arpack_eigsolver<Scalar,form>::class_arpack_eigsolver() {
    RitzToString = {
            {eigsolver_properties::Ritz::LA, "C"},
            {eigsolver_properties::Ritz::SA, "SA"},
            {eigsolver_properties::Ritz::LM, "LM"},
            {eigsolver_properties::Ritz::SM, "SM"},
            {eigsolver_properties::Ritz::LR, "LR"},
            {eigsolver_properties::Ritz::SR, "SR"},
            {eigsolver_properties::Ritz::LI, "LI"},
            {eigsolver_properties::Ritz::SI, "SI"},
            {eigsolver_properties::Ritz::BE, "BE"}
    };
    t_sol.set_properties(profile_eigsolver, 10,"Time iterating  ");
    t_get.set_properties(profile_eigsolver, 10,"Time getting sol");
    t_sub.set_properties(profile_eigsolver, 10,"Time subtracting");
    t_all.set_properties(profile_eigsolver, 10,"Time doing all  ");

}

template<typename Scalar, Form form>
class_arpack_eigsolver<Scalar,form>::class_arpack_eigsolver(
                       const double eigThreshold_,
                       const int eigMaxIter_,
                       const int eigMaxNcv_,
                       const bool compute_eigvecs_,
                       const bool remove_phase_)
            : class_arpack_eigsolver()
{
    setThreshold(eigThreshold_);
    setMaxIter(eigMaxIter_);
    setMaxNcv(eigMaxNcv_);
    compute_eigvecs  = compute_eigvecs_;
    remove_phase = remove_phase_;
}



template<typename Scalar, Form form>
const std::vector<Scalar> & class_arpack_eigsolver<Scalar,form>::ref_eigvecs() const{
    if(eigvecs_found) {
        return eigvecs;
    }else{
        std::cerr << "Eigenvectors haven't been computed yet. Exiting " << std::endl;
        exit(1);
    }
}

template<typename Scalar, Form form>
const std::vector<Scalar> & class_arpack_eigsolver<Scalar,form>::ref_eigvals() const {
    if(eigvals_found) {
        return eigvals;
    }else{
        std::cerr << "Eigenvalues haven't been computed yet. Exiting " << std::endl;
        exit(1);
    }
}


template<typename Scalar, Form form>
const std::vector<Scalar>  class_arpack_eigsolver<Scalar,form>::get_eigvecs() const{
    if(eigvecs_found) {
        return eigvecs;
    }else{
        std::cerr << "Eigenvectors haven't been computed yet. Exiting " << std::endl;
        exit(1);
    }
}

template<typename Scalar, Form form>
const std::vector<Scalar>  class_arpack_eigsolver<Scalar,form>::get_eigvals() const {
    if(eigvals_found) {
        return eigvals;
    }else{
        std::cerr << "Eigenvalues haven't been computed yet. Exiting " << std::endl;
        exit(1);
    }
}


template<typename Scalar,Form form>
void class_arpack_eigsolver<Scalar,form>::shift_invert_eigvals(Scalar sigma) {
    std::transform(eigvals.begin(), eigvals.end(), eigvals.begin(),
                   [sigma](Scalar num) -> Scalar
                   { return 1.0/(num - sigma); });
}

template<typename Scalar, Form form>
void class_arpack_eigsolver<Scalar,form>::subtract_phase() {
    if constexpr (std::is_same<Scalar, std::complex<double>>::value) {
        using namespace std::complex_literals;
        for (int i = 0; i < nev_found; i++) {
            auto begin = eigvecs.begin() + i * rows;
            auto end = begin + rows;
            Scalar inv_phase = -1.0i * std::arg(eigvecs[i * rows]);
            Scalar exp_inv_phase = std::exp(inv_phase);
            std::transform(begin, end, begin,
                           [exp_inv_phase](std::complex<double> num) -> std::complex<double>
                           { return (num * exp_inv_phase); });
        }
    }
}


template<typename Scalar,Form form>
void class_arpack_eigsolver<Scalar,form>::eig(const Scalar *matrix_data,
                                              const int n,
                                              const int nev,
                                              const int ncv,
                                              const Ritz ritz,
                                              const Side side,
                                              const bool compute_eigvecs_,
                                              const bool remove_phase_,
                                              Scalar *residual_
) {
    compute_eigvecs = compute_eigvecs_;
    remove_phase    = remove_phase_;
    residual        = residual_;

    eigvecs_found = false;
    eigvals_found = false;
    eigvals.clear();
    eigvecs.clear();
    int nev_internal = std::min(n/2,nev);
    int ncv_internal = std::max(ncv, 2+nev);
    ncv_internal = std::min(ncv_internal, n);
    assert(ncv_internal >= nev + 2 and ncv_internal <= n);
    assert(nev_internal >= 1 and nev_internal <= n);

    // Stupidly enough, the only difference between the github and "apt" versions of arpack++, is that the apt version only accepts char*
    RitzToString.at(ritz).copy(ritz_char, 2);
    if constexpr(std::is_same_v<Scalar, double> and form == Form::GENERAL) {
        int nev_temp = nev_internal == 1 ? 2 : nev_internal;
        DenseMatrixProduct<Scalar, form> matrix(n, matrix_data, side);
        ARNonSymStdEig<Scalar, DenseMatrixProduct<Scalar, form>> eig(n, nev_temp, &matrix,
                                                                     &DenseMatrixProduct<Scalar, form>::MultMv,
                                                                     ritz_char, ncv_internal, eigThreshold,
                                                                     eigMaxIter, residual);
        counter = matrix.counter;
        find_solution(eig, nev, compute_eigvecs);
    }

    if constexpr(std::is_same_v<Scalar, double> and form == Form::SYMMETRIC) {
        DenseMatrixProduct<Scalar, form> matrix(n, matrix_data, side);
        ARSymStdEig<Scalar, DenseMatrixProduct<Scalar, form>> eig(n, nev_internal, &matrix,
                                                                  &DenseMatrixProduct<Scalar, form>::MultMv,
                                                                  ritz_char, ncv_internal, eigThreshold, eigMaxIter,
                                                                  residual);
        counter = matrix.counter;
        find_solution(eig, nev, compute_eigvecs);
    }

    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        DenseMatrixProduct<Scalar, form> matrix(n, matrix_data, side);
        ARCompStdEig<double, DenseMatrixProduct<Scalar, form>> eig(n, nev_internal, &matrix,
                                                                   &DenseMatrixProduct<Scalar, form>::MultMv,
                                                                   ritz_char, ncv_internal, eigThreshold, eigMaxIter,
                                                                   residual);
        counter = matrix.counter;
        find_solution(eig, nev_internal, compute_eigvecs);
    }
    if (remove_phase) {
        subtract_phase();
    }
}


template<typename Scalar,Form form>
void class_arpack_eigsolver<Scalar,form>::eig_shift_invert(
                                              Scalar *matrix_data,
                                              const int n,
                                              const int nev,
                                              const int ncv,
                                              const Scalar shift,
                                              const Ritz ritz,
                                              const bool compute_eigvecs_,
                                              const bool remove_phase_,
                                              Scalar *residual_
) {
    compute_eigvecs = compute_eigvecs_;
    remove_phase    = remove_phase_;
    residual        = residual_;
    Scalar sigma = shift;
    eigvecs_found = false;
    eigvals_found = false;
    eigvals.clear();
    eigvecs.clear();

    int nev_internal = std::min(n/2,nev);
    int ncv_internal = std::max(ncv, 2+nev);
    ncv_internal = std::min(ncv_internal, n);
    assert(ncv_internal >= nev + 2 and ncv_internal <= n);
    assert(nev_internal >= 1 and nev_internal <= n);

    // Stupidly enough, the only difference between the github and "apt" versions of arpack++, is that the apt version only accepts char*
    RitzToString.at(ritz).copy(ritz_char, 2);
    if constexpr(std::is_same_v<Scalar, double> and form == Form::GENERAL) {
        int nev_temp = nev_internal == 1 ? 2 : nev_internal;
        ARdsNonSymMatrix<Scalar,Scalar> matrix(n, matrix_data);
        ARluNonSymStdEig<Scalar> eig(nev_temp, matrix, sigma,ritz_char, ncv_internal, eigThreshold, eigMaxIter, residual,true);
//        counter = matrix.counter;
        find_solution(eig, nev, compute_eigvecs);
    }

    if constexpr(std::is_same_v<Scalar, double> and form == Form::SYMMETRIC) {
        ARdsSymMatrix<Scalar> matrix(n, matrix_data);
        ARluSymStdEig<Scalar> eig(nev_internal, matrix,sigma, ritz_char, ncv_internal, eigThreshold, eigMaxIter, residual,true);
//        counter = matrix.counter;
        find_solution(eig, nev, compute_eigvecs);
    }

    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        ARdsNonSymMatrix<std::complex<double>,double> matrix(n, matrix_data);
        ARluCompStdEig<double> eig(nev_internal, matrix, sigma,ritz_char, ncv_internal, eigThreshold, eigMaxIter, residual,true);

        find_solution(eig, nev_internal, compute_eigvecs);
//        counter = matrix.counter;
    }

    if (remove_phase) {
        subtract_phase();
    }
}




template<typename Scalar, Form form>
void class_arpack_eigsolver<Scalar,form>::optimize_mps(
        const Scalar *Lblock,        /*!< The left block tensor.  */
        const Scalar *Rblock,        /*!< The right block tensor.  */
        const Scalar *HA,            /*!< The left Hamiltonian MPO's  */
        const Scalar *HB,            /*!< The right Hamiltonian MPO's */
        const std::array<long,4> shape_theta4,         /*!< An array containing the shapes of theta  */
        const std::array<long,4> shape_mpo4 ,           /*!< An array containing the shapes of the MPO  */
        const int nev,
        const int ncv,
        const Ritz ritz,
        const bool remove_phase_,
        Scalar *residual_)
{
    t_all.tic();
    eigvals.clear();
    eigvecs.clear();
    eigvecs_found = false;
    eigvals_found = false;

    remove_phase = remove_phase_;
    residual     = residual_;

    // Stupidly enough, the only difference between the github and "apt" versions of arpack++, is that the apt version only accepts char*
    RitzToString.at(ritz).copy(ritz_char, 2);

//    Lblock->block.data_struct(), Rblock->block.data_struct(), HA->MPO.data_struct(), HB->MPO.data_struct()


    DenseHamiltonianProduct<Scalar>  hamiltonianProduct(Lblock, Rblock,HA, HB, shape_theta4,shape_mpo4);
    int dim  = hamiltonianProduct.cols();
    int ncv_internal = std::max(ncv, 2+nev);
    ncv_internal = std::min(ncv_internal, dim);
    assert(ncv_internal >= nev + 2 and ncv_internal <= dim);

    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        ARCompStdEig<double, DenseHamiltonianProduct<Scalar>> eigsolver(dim, nev, &hamiltonianProduct,
                                                                        &DenseHamiltonianProduct<Scalar>::MultMv, ritz_char,
                                                                        ncv_internal, eigThreshold, eigMaxIter, residual);
        t_sol.tic();
        eigsolver.FindEigenvectors();
        t_sol.toc();
        t_get.tic();
        eigvals_found = eigsolver.EigenvaluesFound();
        eigvecs_found = eigsolver.EigenvectorsFound();
        n         = eigsolver.GetN();
        rows      = eigsolver.GetN();
        cols      = eigsolver.GetNev();
        nev_found = eigsolver.GetNev();
        iter      = eigsolver.GetIter();
        counter   = hamiltonianProduct.counter;
        eigvals   = std::vector<Scalar>(eigsolver.RawEigenvalues() , eigsolver.RawEigenvalues()  + cols);
        eigvecs   = std::vector<Scalar>(eigsolver.RawEigenvector(0), eigsolver.RawEigenvectors() + rows*cols);
        t_get.toc();
    }

    if constexpr(std::is_same_v<Scalar, double>) {
        ARNonSymStdEig <double, DenseHamiltonianProduct<Scalar>> eigsolver(dim, nev, &hamiltonianProduct,
                                                                           &DenseHamiltonianProduct<Scalar>::MultMv, ritz_char,
                                                                           ncv_internal, eigThreshold, eigMaxIter, residual);
        eigsolver.FindEigenvectors();
        eigvals_found   = eigsolver.EigenvaluesFound();
        eigvecs_found   = eigsolver.EigenvectorsFound();
        n               = eigsolver.GetN();
        rows            = eigsolver.GetN();
        cols            = eigsolver.GetNev();
        nev_found       = eigsolver.GetNev();
        iter            = eigsolver.GetIter();
        counter         = hamiltonianProduct.counter;
        eigvals         = std::vector<Scalar>(eigsolver.RawEigenvalues() , eigsolver.RawEigenvalues()  + cols);
        eigvecs         = std::vector<Scalar>(eigsolver.RawEigenvector(0), eigsolver.RawEigenvectors() + rows*cols);
    }

    if (remove_phase_){
        t_sub.tic();
        subtract_phase();
        t_sub.toc();
    }
    t_all.toc();
    hamiltonianProduct.t_mul.print_time();
    t_sol.print_time();
    t_get.print_time();
    t_sub.print_time();
    t_all.print_time();

}




// Explicit instantiations


template class class_arpack_eigsolver<std::complex<double>, Form::GENERAL> ;
template class class_arpack_eigsolver<std::complex<double>, Form::SYMMETRIC> ;

template class class_arpack_eigsolver<double, Form::GENERAL> ;
template class class_arpack_eigsolver<double, Form::SYMMETRIC> ;
