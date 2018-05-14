//
// Created by david on 2018-05-06.
//

#include "class_arpack_eigsolver.h"
#include <assert.h>
#include <algorithm>
#include <arpackpp/ardsnsym.h>
#include <arpackpp/ardscomp.h>
#include <arpackpp/ardgcomp.h>
#include <arpackpp/ardssym.h>
#include <general/class_arpack_custom_products.h>
#include <sim_parameters/nmspc_sim_settings.h>


template<typename Scalar, Form form>
class_arpack_eigsolver<Scalar,form>::class_arpack_eigsolver() {
    RitzToString = {
            {arpackpp::Ritz::LA, "LA"},
            {arpackpp::Ritz::SA, "SA"},
            {arpackpp::Ritz::LM, "LM"},
            {arpackpp::Ritz::SM, "SM"},
            {arpackpp::Ritz::LR, "LR"},
            {arpackpp::Ritz::SR, "SR"},
            {arpackpp::Ritz::LI, "LI"},
            {arpackpp::Ritz::SI, "SI"},
            {arpackpp::Ritz::BE, "BE"}
    };
}
template<typename Scalar, Form form>
class_arpack_eigsolver<Scalar, form>::class_arpack_eigsolver(Scalar *matrix_data,
                                                                     Ritz ritz,
                                                                     Side side,
                                                                     const int n,
                                                                     const int nev,
                                                                     const int ncv,
                                                                     bool getvecs,
                                                                     bool dephase,
                                                                     Scalar *residp)
: class_arpack_eigsolver()
{
    eig(matrix_data,ritz,side,n,nev,ncv,getvecs,dephase,residp);
}
template<typename Scalar, Form form>
class_arpack_eigsolver<Scalar, form>::class_arpack_eigsolver(
        const Scalar *Lblock,                           /*!< The left block tensor.  */
        const Scalar *Rblock,                           /*!< The right block tensor.  */
        const Scalar *HA,                               /*!< The left Hamiltonian MPO's  */
        const Scalar *HB,                               /*!< The right Hamiltonian MPO's */
        const std::array<long,4> shape_theta4,          /*!< An array containing the shapes of theta  */
        const std::array<long,4> shape_mpo4 ,           /*!< An array containing the shapes of the MPO  */
        Ritz ritz,
        int nev,
        int ncv,
        bool bool_dephase,
        Scalar *resid)
: class_arpack_eigsolver()

{
 optimize_mps(Lblock,Rblock, HA,HB, shape_theta4, shape_mpo4, ritz, nev, ncv, bool_dephase, resid);
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





template<typename Scalar, Form form>
void class_arpack_eigsolver<Scalar,form>::subtract_phase() {
    if constexpr (std::is_same_v<Scalar, std::complex<double>>) {
        using namespace std::complex_literals;
        for (int i = 0; i < nev_found; i++) {
            auto begin = eigvecs.begin() + i * rows;
            auto end = begin + rows;
            Scalar inv_phase = -1.0i * std::arg(eigvecs[i * rows]);
            Scalar exp_inv_phase = std::exp(inv_phase);
//            std::transform(begin, end, begin,
//                           std::bind1st(std::multiplies<Scalar>(), exp_inv_phase));
            std::transform(begin, end, begin,
                           [exp_inv_phase](std::complex<double> num) -> std::complex<double>
                                   { return (num * exp_inv_phase).real(); });
        }
    }
}


template<typename Scalar,Form form>
void class_arpack_eigsolver<Scalar,form>::eig(Scalar *matrix_data,
                                          Ritz ritz,
                                          Side side,
                                          const int n,
                                          const int nev,
                                          const int ncv,
                                          bool bool_find_eigvecs,
                                          bool bool_dephase,
                                          Scalar *residp
) {

    eigvecs_found = false;
    eigvals_found = false;
    eigvals.clear();
    eigvecs.clear();

    if constexpr(std::is_same_v<Scalar, double> and form == Form::GENERAL) {
        int nev_temp = nev == 1 ? 2 : nev;
        DenseMatrixProduct<Scalar, form> matrix(n, matrix_data, side);
        ARNonSymStdEig<Scalar, DenseMatrixProduct<Scalar, form>> eig(n, nev_temp, &matrix,
                                                                     &DenseMatrixProduct<Scalar, form>::MultMv,
                                                                     RitzToString.at(ritz), ncv, eigThreshold,
                                                                     eigMaxIter, residp);
        counter = matrix.counter;
        find_solution(eig, nev, bool_find_eigvecs);
    }

    if constexpr(std::is_same_v<Scalar, double> and form == Form::SYMMETRIC) {
        DenseMatrixProduct<Scalar, form> matrix(n, matrix_data, side);
        ARSymStdEig<Scalar, DenseMatrixProduct<Scalar, form>> eig(n, nev, &matrix,
                                                                  &DenseMatrixProduct<Scalar, form>::MultMv,
                                                                  RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter,
                                                                  residp);
        counter = matrix.counter;
        find_solution(eig, nev, bool_find_eigvecs);
    }

    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        DenseMatrixProduct<Scalar, form> matrix(n, matrix_data, side);
        ARCompStdEig<double, DenseMatrixProduct<Scalar, form>> eig(n, nev, &matrix,
                                                                   &DenseMatrixProduct<Scalar, form>::MultMv,
                                                                   RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter,
                                                                   residp);
        counter = matrix.counter;
        find_solution(eig, nev, bool_find_eigvecs);
    }
    if (bool_dephase) {
        subtract_phase();
    }
}



template<typename Scalar, Form form>
const std::pair<const std::vector<Scalar>&,const std::vector<Scalar>&>
class_arpack_eigsolver<Scalar,form>::eig_ref_vec_val(Scalar *matrix_data,
                                                      Ritz ritz,
                                                      Side side,
                                                      const int n,
                                                      const int nev,
                                                      const int ncv,
                                                      bool bool_dephase,
                                                      Scalar *residp
)
{
    eig(matrix_data,ritz,side, n, nev, ncv,true,bool_dephase,residp);
    return make_pair(ref_eigvecs(), ref_eigvals());
};

template<typename Scalar, Form form>
const std::pair<const std::vector<Scalar>,const std::vector<Scalar>>
class_arpack_eigsolver<Scalar,form>::eig_get_vec_val(Scalar *matrix_data,
                                                      Ritz ritz,
                                                      Side side,
                                                      const int n,
                                                      const int nev,
                                                      const int ncv,
                                                      bool bool_dephase,
                                                      Scalar *residp
)
{
    eig(matrix_data,ritz,side, n, nev, ncv,true,bool_dephase,residp);
    return make_pair(get_eigvecs(), get_eigvals());
};


template<typename Scalar, Form form>
void class_arpack_eigsolver<Scalar,form>::optimize_mps(
        const Scalar *Lblock,        /*!< The left block tensor.  */
        const Scalar *Rblock,        /*!< The right block tensor.  */
        const Scalar *HA,            /*!< The left Hamiltonian MPO's  */
        const Scalar *HB,            /*!< The right Hamiltonian MPO's */
        const std::array<long,4> shape_theta4,         /*!< An array containing the shapes of theta  */
        const std::array<long,4> shape_mpo4 ,           /*!< An array containing the shapes of the MPO  */
        Ritz ritz,
        int nev,
        int ncv,
        bool bool_dephase,
        Scalar *resid)
{
    eigvals.clear();
    eigvecs.clear();
    eigvecs_found = false;
    eigvals_found = false;
    DenseHamiltonianProduct<Scalar>  hamiltonianProduct(Lblock, Rblock, HA, HB, shape_theta4,shape_mpo4);
    int dim  = hamiltonianProduct.cols();
//    int size = hamiltonianProduct.cols() * hamiltonianProduct.rows();
    int ncv_internal = std::max(dim/2, 4*nev);
    ncv_internal = std::min(ncv, ncv_internal);

    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        ARCompStdEig<double, DenseHamiltonianProduct<Scalar>> eigsolver(dim, nev, &hamiltonianProduct,
                                                                        &DenseHamiltonianProduct<Scalar>::MultMv, RitzToString.at(ritz),
                                                                        ncv_internal, eigThreshold, eigMaxIter, resid);

        eigsolver.FindEigenvectors();
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

    }

    if constexpr(std::is_same_v<Scalar, double>) {
        ARNonSymStdEig <double, DenseHamiltonianProduct<Scalar>> eigsolver(dim, nev, &hamiltonianProduct,
                                                                        &DenseHamiltonianProduct<Scalar>::MultMv, RitzToString.at(ritz),
                                                                        ncv_internal, eigThreshold, eigMaxIter, resid);

        eigsolver.FindEigenvectors();
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

    }

   if (bool_dephase){
        subtract_phase();
    }
}




// Explicit instantiations


template class class_arpack_eigsolver<std::complex<double>, Form::GENERAL> ;
template class class_arpack_eigsolver<std::complex<double>, Form::SYMMETRIC> ;

template class class_arpack_eigsolver<double, Form::GENERAL> ;
template class class_arpack_eigsolver<double, Form::SYMMETRIC> ;

