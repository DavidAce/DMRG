//
// Created by david on 2018-10-29.
//

#ifndef EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H
#define EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H

#include "arpack_extra/matrix_product_dense.h"
#include "arpack_extra/matrix_product_sparse.h"
#include "arpack_extra/arpackpp_solver.h"
#include "general/nmspc_eigutils.h"
#include "arpack_extra/matrix_recast.h"


class class_eigsolver {
public:

    eigutils::eigConfig solverConf;
    eigutils::eigSolution   solution;

    void conf_init(const int L,
                   const int nev,
                   const int ncv,
                   const std::complex<double> sigma            = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                   const eigutils::eigSetting::Type type        = eigutils::eigSetting::Type::REAL,
                   const eigutils::eigSetting::Form form        = eigutils::eigSetting::Form::SYMMETRIC,
                   const eigutils::eigSetting::Ritz ritz        = eigutils::eigSetting::Ritz::LM,
                   const eigutils::eigSetting::Side side        = eigutils::eigSetting::Side::R,
                   const eigutils::eigSetting::Storage storage  = eigutils::eigSetting::Storage::DENSE,
                   const bool compute_eigvecs_           = false,
                   const bool remove_phase_              = false
                   );


    template<typename Scalar>
    void eigs_auto(const Scalar *matrix_data,
                   const int L,
                   const int nev,
                   const bool compute_eigvecs_           = false,
                   const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
                   const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                   const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                   const bool remove_phase_              = false);


    template<eigutils::eigSetting::Storage storage,typename Scalar>
    void eigs (const  Scalar *matrix,
               const int L,
               const int nev,
               const int ncv,
               const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
               const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
               const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
               const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
               const bool compute_eigvecs_           = false,
               const bool remove_phase_              = false,
               Scalar *residual_                     = nullptr);


    template<typename Scalar>
    void eigs_dense(const Scalar *matrix,
                   const int L,
                   const int nev,
                   const int ncv,
                   const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                   const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
                   const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
                   const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                   const bool compute_eigvecs_           = false,
                   const bool remove_phase_              = false,
                   Scalar *residual_                     = nullptr);


    template<typename Scalar>
    void eigs_dense(DenseMatrixProduct<Scalar> &matrix,
                    const int nev,
                    const int ncv,
                    const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                    const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
                    const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
                    const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                    const bool compute_eigvecs_           = false,
                    const bool remove_phase_              = false,
                    Scalar *residual_                     = nullptr);



    template<typename Scalar>
    void eigs_sparse(const Scalar *matrix,
                     const int L,
                     const int nev,
                     const int ncv,
                     const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                     const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
                     const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
                     const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                     const bool compute_eigvecs_           = false,
                     const bool remove_phase_              = false,
                     Scalar *residual_                     = nullptr);


    template<typename Scalar>
    void eigs_sparse(SparseMatrixProduct<Scalar> &matrix,
                    const int nev,
                    const int ncv,
                    const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                    const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
                    const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::SR,
                    const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                    const bool compute_eigvecs_           = false,
                    const bool remove_phase_              = false,
                    Scalar *residual_                     = nullptr);




};


// Definitions



template<typename Scalar>
void class_eigsolver::eigs_auto   (const Scalar *matrix_data,
                                   const int L,
                                   const int nev,
                                   const bool compute_eigvecs_           ,
                                   const eigutils::eigSetting::Ritz ritz ,
                                   const std::complex<double> sigma      ,
                                   const eigutils::eigSetting::Side side ,
                                   const bool remove_phase_              )
{
    using namespace eigutils::eigSetting;
    matrix_recast<Scalar> matRecast(matrix_data,L);
    bool is_sparse    = matRecast.is_sparse();
    bool is_real      = matRecast.is_real();
    bool is_symmetric = matRecast.is_symmetric();

    Form form        = is_symmetric ? Form::SYMMETRIC : Form::NONSYMMETRIC;
    Type type        = is_real      ? Type::REAL      : Type ::CPLX;
    Storage storage  = is_sparse    ? Storage::SPARSE : Storage::DENSE;

    conf_init(L, nev, -1, sigma,type, form,ritz,side,storage,compute_eigvecs_,remove_phase_);

    if(is_real) {
        if(is_sparse) {
            auto matrix = matRecast.get_as_real_sparse();
            arpackpp_solver<SparseMatrixProduct<double>> solver(matrix, solverConf, solution);
            solver.eigs();
        }else {
            auto matrix = matRecast.get_as_real_dense();
            arpackpp_solver<DenseMatrixProduct<double>> solver(matrix, solverConf, solution);
            solver.eigs();
        }
    }else {
        if(is_sparse) {
            auto matrix = matRecast.get_as_cplx_sparse();
            arpackpp_solver<SparseMatrixProduct<std::complex<double>>> solver(matrix, solverConf, solution);
            solver.eigs();
        }else {
            auto matrix = matRecast.get_as_cplx_dense();
            arpackpp_solver<DenseMatrixProduct<std::complex<double>>> solver(matrix, solverConf, solution);
            solver.eigs();
        }
    }
}





template<eigutils::eigSetting::Storage storage,typename Scalar>
void class_eigsolver::eigs (const Scalar *matrix,
                            const int L,
                            const int nev,
                            const int ncv,
                            const std::complex<double> sigma,
                            const eigutils::eigSetting::Form form,
                            const eigutils::eigSetting::Ritz ritz,
                            const eigutils::eigSetting::Side side,
                            const bool compute_eigvecs_,
                            const bool remove_phase_,
                            Scalar *residual_)
{
    using namespace eigutils::eigSetting;
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    if constexpr(storage == Storage::DENSE){
        auto matrix_dense = DenseMatrixProduct<Scalar> (matrix,L);
        arpackpp_solver<DenseMatrixProduct<Scalar>> solver(matrix_dense, solverConf, solution,residual_);
        solver.eigs();
    }else{
        auto matrix_sparse = SparseMatrixProduct<Scalar> (matrix,L);
        arpackpp_solver<SparseMatrixProduct<Scalar>> solver(matrix_sparse, solverConf, solution,residual_);
        solver.eigs();
    }
}



template<typename Scalar>
void class_eigsolver::eigs_dense   (const Scalar *matrix,
                                   const int L,
                                   const int nev,
                                   const int ncv,
                                   const std::complex<double> sigma,
                                   const eigutils::eigSetting::Form form,
                                   const eigutils::eigSetting::Ritz ritz,
                                   const eigutils::eigSetting::Side side,
                                   const bool compute_eigvecs_,
                                   const bool remove_phase_,
                                   Scalar *residual_)
{
    using namespace eigutils::eigSetting;
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    Storage storage = Storage::DENSE;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    auto matrix_dense = DenseMatrixProduct<Scalar> (matrix,L);
    arpackpp_solver<DenseMatrixProduct<Scalar>> solver(matrix_dense, solverConf, solution,residual_);
    solver.eigs();
}


template<typename Scalar>
void class_eigsolver::eigs_dense   (DenseMatrixProduct<Scalar> &matrix_dense,
                                    const int nev,
                                    const int ncv,
                                    const std::complex<double> sigma,
                                    const eigutils::eigSetting::Form form,
                                    const eigutils::eigSetting::Ritz ritz,
                                    const eigutils::eigSetting::Side side,
                                    const bool compute_eigvecs_,
                                    const bool remove_phase_,
                                    Scalar *residual_)
{
    using namespace eigutils::eigSetting;
    int L = matrix_dense.rows();
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    Storage storage = Storage::DENSE;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    arpackpp_solver<DenseMatrixProduct<Scalar>> solver(matrix_dense, solverConf, solution,residual_);
    solver.eigs();
}



template<typename Scalar>
void class_eigsolver::eigs_sparse   (const Scalar *matrix,
                                    const int L,
                                    const int nev,
                                    const int ncv,
                                    const std::complex<double> sigma,
                                    const eigutils::eigSetting::Form form,
                                    const eigutils::eigSetting::Ritz ritz,
                                    const eigutils::eigSetting::Side side,
                                    const bool compute_eigvecs_,
                                    const bool remove_phase_,
                                    Scalar *residual_)
{
    using namespace eigutils::eigSetting;
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    Storage storage = Storage::SPARSE;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    auto matrix_sparse = SparseMatrixProduct<Scalar> (matrix,L);
    arpackpp_solver<SparseMatrixProduct<Scalar>> solver(matrix_sparse, solverConf, solution,residual_);
    solver.eigs();
}




template<typename Scalar>
void class_eigsolver::eigs_sparse   (SparseMatrixProduct<Scalar> &matrix_sparse,
                                     const int nev,
                                     const int ncv,
                                     const std::complex<double> sigma,
                                     const eigutils::eigSetting::Form form,
                                     const eigutils::eigSetting::Ritz ritz,
                                     const eigutils::eigSetting::Side side,
                                     const bool compute_eigvecs_,
                                     const bool remove_phase_,
                                     Scalar *residual_)
{
    using namespace eigutils::eigSetting;
    int L = matrix_sparse.rows();
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    Storage storage = Storage::SPARSE;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    arpackpp_solver<SparseMatrixProduct<Scalar>> solver(matrix_sparse, solverConf, solution,residual_);
    solver.eigs();
}


#endif //EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H
