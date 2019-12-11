//
// Created by david on 2018-10-30.
//

#if  defined(MKL_AVAILABLE)
#include <mkl_lapacke.h>
#elif defined(OpenBLAS_AVAILABLE)
#include <lapacke.h>
#endif

#include "arpackpp_solver.h"

#if __has_include(<arpackpp/arcomp.h>)
#include <arpackpp/arcomp.h>
#include <arpackpp/ardscomp.h>
#include <arpackpp/ardnsmat.h>
#include <arpackpp/arsnsym.h>
#include <arpackpp/arssym.h>
#elif __has_include(<arpack++/arcomp.h>)
#include <arpack++/arcomp.h>
#include <arpack++/ardscomp.h>
#include <arpack++/ardnsmat.h>
#include <arpack++/arsnsym.h>
#include <arpack++/arssym.h>
#else
#error Could not include arpack headers correctly
#endif

#include "matrix_product_dense.h"
#include "matrix_product_sparse.h"
#include "matrix_product_stl.h"
#include "matrix_product_hamiltonian.h"
#include <general/nmspc_type_check.h>


namespace tc = TypeCheck;
using namespace eigutils::eigSetting;


template<typename MatrixType>
arpackpp_solver<MatrixType>::arpackpp_solver(
        MatrixType               &matrix_,
        eigutils::eigConfig      &solverConf_,
        eigutils::eigSolution    &solution_,
        Scalar                   *residual_
        )
        :
        matrix(matrix_),
        solverConf(solverConf_),
        solution  (solution_),
        residual  (residual_)
{

    t_sol.set_properties(profile_arpack, 10,"Time iterating  ");
    t_get.set_properties(profile_arpack, 10,"Time getting sol");
    t_sub.set_properties(profile_arpack, 10,"Time subtracting");
    t_all.set_properties(profile_arpack, 10,"Time doing all  ");
}




template<typename MatrixType>
void arpackpp_solver<MatrixType>::eigs() {
    solution.reset();
    nev_internal = std::min(matrix.rows()/2,solverConf.eigMaxNev);
    ncv_internal = std::max(solverConf.eigMaxNcv, 2+solverConf.eigMaxNev);
    ncv_internal = std::min(ncv_internal, matrix.rows());
    assert(ncv_internal >= solverConf.eigMaxNev + 2 and ncv_internal <= matrix.rows());
    assert(nev_internal >= 1 and nev_internal <= matrix.rows() / 2);

    // Stupidly enough, the only difference between the github and "apt" versions of arpack++,
    // is that the apt version only accepts char*, whereas the github one accepts string and char*.
    // For this reason we have to convert the ritz to a format that both can take.
    solverConf.writeRitzChar();
    matrix.set_mode(solverConf.form);
    matrix.set_side(solverConf.side);

    // Calculate shift-inverse mat-vec mult operator by LU decomposition
    if constexpr(MatrixType::can_shift){
        if(solverConf.shift == eigutils::eigSetting::Shift::ON) {
            matrix.set_shift(solverConf.sigma);
            matrix.FactorOP();
        }
    }
    if(not solverConf.confOK) throw std::runtime_error("solverConf isn't ready!");
    assert(solverConf.confOK and "solverConf isn't ready!");
    // Dispatch to symmetric or nonsymmetric. If complex, there's only a nonsymmetric option available.
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        this->eigs_comp();
    }else{
        if(solverConf.form == Form::SYMMETRIC){this->eigs_sym();}else {this->eigs_nsym();}
    }

}



template<typename MatrixType>
void arpackpp_solver<MatrixType>::eigs_sym() {
    if constexpr(std::is_same<Scalar,double>::value) {
        assert(solverConf.form       == Form::SYMMETRIC and "ERROR: solverConf not SYMMETRIC");
        assert(matrix.get_form()     == Form::SYMMETRIC and "ERROR: matrix not SYMMETRIC");
        ARSymStdEig<double, MatrixType> solver(
                matrix.rows(),
                nev_internal,
                &matrix,
                &MatrixType::MultAx,
                solverConf.ritz_char,
                ncv_internal,
                solverConf.eigThreshold,
                solverConf.eigMaxIter,
                residual);
        if constexpr(MatrixType::can_shift){
            switch (solverConf.shift) {
                case Shift::OFF :
                    break;
                case Shift::ON :
                    solver.SetShiftInvertMode(std::real(solverConf.sigma), &matrix, &MatrixType::MultOPv);
                    break;
            }
        }


        this->find_solution(solver, solverConf.eigMaxNev);
        this->copy_solution<Type::REAL, Form::SYMMETRIC>(solver);
    }else{
        eigutils::eigLogger::log->critical("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()) );
        throw std::runtime_error("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}


template<typename MatrixType>
void arpackpp_solver<MatrixType>::eigs_nsym() {
    if constexpr(std::is_same<Scalar, double>::value) {
        assert(solverConf.form == Form::NONSYMMETRIC and "ERROR: solverConf not NONSYMMETRIC");
        assert(matrix.get_form() == Form::NONSYMMETRIC and "ERROR: matrix not NONSYMMETRIC");
        if (nev_internal == 1) { nev_internal++; }

        ARNonSymStdEig<double, MatrixType> solver(
                matrix.rows(),
                nev_internal,
                &matrix,
                &MatrixType::MultAx,
                solverConf.ritz_char,
                ncv_internal,
                solverConf.eigThreshold,
                solverConf.eigMaxIter,
                residual);

        if constexpr(MatrixType::can_shift){
            switch (solverConf.shift) {
                case Shift::OFF :
                    break;
                case Shift::ON :
                    solver.SetShiftInvertMode(std::real(solverConf.sigma), &matrix, &MatrixType::MultOPv);
                    break;
            }
        }

        this->find_solution(solver, solverConf.eigMaxNev);
        if (matrix.get_side() == Side::R){
            this->copy_solution<Type::REAL, Form::NONSYMMETRIC, Side::R>(solver);
        }else{
            this->copy_solution<Type::REAL, Form::NONSYMMETRIC, Side::L>(solver);
        }


    }else{
        eigutils::eigLogger::log->critical("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()) );
        throw std::runtime_error("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}





template<typename MatrixType>
void arpackpp_solver<MatrixType>::eigs_comp() {
    if constexpr(std::is_same<Scalar, std::complex<double>>::value){
        ARCompStdEig<double, MatrixType> solver(
                matrix.rows(),
                nev_internal,
                &matrix,
                &MatrixType::MultAx,
                solverConf.ritz_char,
                ncv_internal,
                solverConf.eigThreshold,
                solverConf.eigMaxIter,
                residual);

        if constexpr(MatrixType::can_shift){
            switch (solverConf.shift) {
                case Shift::OFF :
                    break;
                case Shift::ON :
                    solver.SetShiftInvertMode(std::real(solverConf.sigma), &matrix, &MatrixType::MultOPv);
                    break;
            }
        }
        this->find_solution(solver, solverConf.eigMaxNev);

        if (matrix.get_form() == Form::SYMMETRIC and matrix.get_side() == Side::R ){
            this->copy_solution<Type::CPLX, Form::SYMMETRIC, Side::R>(solver);
        }else if (matrix.get_form() == Form::SYMMETRIC and matrix.get_side() == Side::L ){
            this->copy_solution<Type::CPLX, Form::SYMMETRIC, Side::L>(solver);
        }else if (matrix.get_form() == Form::NONSYMMETRIC and matrix.get_side() == Side::R ){
            this->copy_solution<Type::CPLX, Form::NONSYMMETRIC, Side::R>(solver);
        }else if (matrix.get_form() == Form::NONSYMMETRIC and matrix.get_side() == Side::L ){
            this->copy_solution<Type::CPLX, Form::NONSYMMETRIC, Side::L>(solver);
        }


    }else{
        eigutils::eigLogger::log->critical("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()) );
        throw std::runtime_error("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}


//
//
template <typename MatrixType>
template <typename Derived>
void arpackpp_solver<MatrixType>::find_solution(Derived &solver, int nev) {
    if (solverConf.compute_eigvecs) {
        solver.FindEigenvectors();
        solution.meta.eigvals_found  = solver.EigenvaluesFound();  //BOOL!
        solution.meta.eigvecsR_found = solver.EigenvectorsFound(); //BOOL!
        solution.meta.iter           = solver.GetIter();
        solution.meta.n              = solver.GetN();
        solution.meta.nev            = std::min(nev, solver.GetNev());
        solution.meta.nev_converged  = solver.ConvergedEigenvalues();
        solution.meta.ncv_used       = solver.GetNcv();
        solution.meta.rows           = solver.GetN();
        solution.meta.cols           = solution.meta.nev;
        solution.meta.counter        = matrix.counter;
    }else{
        solver.FindEigenvalues();
        solution.meta.eigvals_found = solver.EigenvaluesFound();
        solution.meta.iter          = solver.GetIter();
        solution.meta.n             = solver.GetN();
        solution.meta.nev           = std::min(nev, solver.GetNev());
        solution.meta.nev_converged = solver.ConvergedEigenvalues();
        solution.meta.ncv_used      = solver.GetNcv();
        solution.meta.rows          = solver.GetN();
        solution.meta.cols          = solution.meta.nev;
        solution.meta.counter       = matrix.counter;
    }
}

//
//
//template <typename MatrixType>
//template <typename Derived>
//void arpackpp_solver<MatrixType>::copy_solution_symm(Derived &solver) {
//    int eigvecsize = solution.meta.n * solution.meta.nev;
//    int eigvalsize = solution.meta.nev;
//    if (solverConf.compute_eigvecs) {
//        if(solverConf.side == Side::R){
//            if constexpr(std::is_same<double,typename MatrixType::Scalar>::value ){
//                solution.eigvecsR_real.resize(eigvecsize);
//                solution.eigvals_real.resize(eigvalsize);
//                std::copy(solver.RawEigenvectors(),solver.RawEigenvectors() + eigvecsize, solution.eigvecsR_real.begin());
//                std::copy(solver.RawEigenvalues() , solver.RawEigenvalues() + eigvalsize, solution.eigvals_real.begin());
//            }else if constexpr(std::is_same<std::complex<double>,typename MatrixType::Scalar>::value ){
//                solution.eigvecsR_cplx.resize(eigvecsize);
//                solution.eigvals_cplx.resize(eigvalsize);
//                std::copy(solver.RawEigenvectors(),solver.RawEigenvectors() + eigvecsize, solution.eigvecsR_cplx.begin());
//                std::copy(solver.RawEigenvalues() ,solver.RawEigenvalues() + eigvalsize, solution.eigvals_cplx.begin());
//                if(solverConf.form == Form::SYMMETRIC){
//                    for (int j = 0; j < eigvalsize; j++) {
//                        solution.eigvals_real.emplace_back(solution.eigvals_cplx[j].real());
//                    }
//                }
//            }
//        }
//        if(solverConf.side == Side::L){
//            if constexpr(std::is_same<double,typename MatrixType::Scalar>::value ){
//                solution.eigvecsL_real.resize(eigvecsize);
//                solution.eigvals_real.resize(eigvalsize);
//                std::copy(solver.RawEigenvectors(),solver.RawEigenvectors() + eigvecsize, solution.eigvecsL_real.begin());
//                std::copy(solver.RawEigenvalues() ,solver.RawEigenvalues() + eigvalsize, solution.eigvals_real.begin());
//
//            }else if constexpr(std::is_same<std::complex<double>,typename MatrixType::Scalar>::value ){
//                solution.eigvecsL_cplx.resize(eigvecsize);
//                solution.eigvals_cplx.resize(eigvalsize);
//                std::copy(solver.RawEigenvectors(),solver.RawEigenvectors() + eigvecsize, solution.eigvecsL_cplx.begin());
//                std::copy(solver.RawEigenvalues() ,solver.RawEigenvalues() + eigvalsize, solution.eigvals_cplx.begin());
//                if(solverConf.form == Form::SYMMETRIC){
//                    for (int j = 0; j < eigvalsize; j++) {
//                        solution.eigvals_real.emplace_back(solution.eigvals_cplx[j].real());
//                    }
//                }
//            }
//
//        }
//    }
//
//}


//
//template <typename MatrixType>
//template <typename Derived>
//void arpackpp_solver<MatrixType>::copy_solution_nsym(Derived &solver) {
//    for (int j = 0; j < solution.meta.cols; j++) {
//        solution.eigvals_cplx.emplace_back(std::complex<double>(solver.EigenvalueReal(j), solver.EigenvalueImag(j)));
//    }
//    if(solverConf.compute_eigvecs){
//        for (int j = 0; j < solution.meta.cols; j++){
//            for (int i = 0; i < solution.meta.rows; i++){
//                solution.eigvecsR_cplx.emplace_back(std::complex<double>(solver.EigenvectorReal(j,i), solver.EigenvectorImag(j,i)));
//            }
//        }
//    }
//}




// Explicit instantiations
//
template class arpackpp_solver<DenseMatrixProduct<double>>;
template class arpackpp_solver<DenseMatrixProduct<std::complex<double>>>;
template class arpackpp_solver<SparseMatrixProduct<double>>;
template class arpackpp_solver<SparseMatrixProduct<std::complex<double>>>;
template class arpackpp_solver<StlMatrixProduct<double>>;
template class arpackpp_solver<StlMatrixProduct<std::complex<double>>>;
template class arpackpp_solver<DenseHamiltonianProduct<double>>;
template class arpackpp_solver<DenseHamiltonianProduct<std::complex<double>>>;
