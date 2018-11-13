//
// Created by david on 2018-05-08.
//

#ifndef MATRIX_PRODUCT_SPARSE_H
#define MATRIX_PRODUCT_SPARSE_H

#ifdef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS_SUSPEND
#undef EIGEN_USE_BLAS
#endif


#include <general/nmspc_eigutils.h>
#include <general/class_tic_toc.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#define profile_matrix_product_sparse 1


template <typename Scalar_>
class SparseMatrixProduct {
public:
    using Scalar          = Scalar_;
    using MatrixType      = Eigen::SparseMatrix<Scalar>;
    using DenseMatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using VectorType      = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    using VectorTypeT     = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;
private:
    MatrixType A_matrix;          // The actual matrix. Given matrices will be copied into this one.
    const int L;                        // The linear matrix dimension
    eigutils::eigSetting::Form form;     // Chooses SYMMETRIC / NONSYMMETRIC mode
    eigutils::eigSetting::Side side;     // Chooses whether to find (R)ight or (L)eft eigenvectors

    // Shift-invert mode stuff
//    MatrixType A_shift;
    Eigen::SparseLU<MatrixType>          lu_sparse;              // Object for sparse LU decomposition used in shift-invert mode
    Eigen::PartialPivLU<DenseMatrixType> lu_dense;               // Object for dense LU decomposition used in shift-invert mode
//    Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Eigen::Index>, Eigen::COLAMDOrdering<Eigen::Index> > lu_sparse;

    double sigmaR = std::numeric_limits<double>::quiet_NaN();   // The real part of the shift
    double sigmaI = std::numeric_limits<double>::quiet_NaN();   // The imag part of the shift
    bool sparseLU      = false;
    bool readyFactorOp = false;                                 // Flag to make sure LU factorization has occurred
    bool readyShift = false;

public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    SparseMatrixProduct(
            const Scalar * A_,
            const int L_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R)
            : A_matrix(Eigen::Map<const DenseMatrixType>(A_,L_,L_).sparseView()), L(L_), form(form_), side(side_)
    {
        A_matrix.makeCompressed();
        init_profiling();
    }

    // Eigen type constructor. Pass any copy-assignable eigen type into an internal Eigen matrix.
    template<typename Derived>
    explicit SparseMatrixProduct(
            const Eigen::EigenBase<Derived> &matrix_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R)
            : A_matrix(matrix_.derived().sparseView()), L(A_matrix.rows()), form(form_), side(side_)
    {
        A_matrix.makeCompressed();
        init_profiling();
    }

    // Functions used in in Arpack++ solver
    int rows() const {return L;};
    int cols() const {return L;};
    void FactorOP();                                      //  Factors (A-sigma*I) into PLU
    void MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr);    //   Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultAx (Scalar* x_in_ptr, Scalar* x_out_ptr);    //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int counter = 0;
    void print()const;
    const auto & get_matrix(){return A_matrix;}
    void set_shift(std::complex<double> sigma_)   {if(readyShift){return;} sigmaR=std::real(sigma_);sigmaI=std::imag(sigma_) ;readyShift = true;}
    void set_shift(double               sigma_)   {if(readyShift){return;} sigmaR=sigma_, sigmaI = 0.0;readyShift = true;}
    void set_shift(double sigmaR_, double sigmaI_){if(readyShift){return;} sigmaR=sigmaR_;sigmaI = sigmaI_ ;readyShift = true;}
    void set_mode(const eigutils::eigSetting::Form form_){form = form_;}
    void set_side(const eigutils::eigSetting::Side side_){side = side_;}
    const MatrixType & get_matrix()const {return A_matrix;}
    const eigutils::eigSetting::Form &get_form()const{return form;}
    const eigutils::eigSetting::Side &get_side()const{return side;}

    // Profiling
    void init_profiling(){
        t_factorOp.set_properties(profile_matrix_product_sparse, 5,"Time FactorOp");
        t_multOpv.set_properties(profile_matrix_product_sparse, 5,"Time MultOpv");
        t_multax.set_properties(profile_matrix_product_sparse, 5,"Time MultAx");
    }
    class_tic_toc t_factorOp;
    class_tic_toc t_multOpv;
    class_tic_toc t_multax;

};




// Function definitions



template<typename Scalar>
void SparseMatrixProduct<Scalar>::print() const {
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}


template<typename Scalar>
void SparseMatrixProduct<Scalar>::FactorOP()

/*  Sparse decomposition
 *  Factors P(A-sigma*I) = LU
 */

{   if(readyFactorOp){return;}
    assert(readyShift and "Shift value sigma has not been set.");
    Scalar sigma;
    t_factorOp.tic();
    MatrixType A_shift;
    if constexpr(std::is_same<Scalar,double>::value)
    {
        sigma = sigmaR;
        A_shift = (A_matrix - sigmaR * Eigen::MatrixXd::Identity(L,L));
    }
    else
    {
        sigma = std::complex<double>(sigmaR,sigmaI);
        A_shift = (A_matrix - sigma * Eigen::MatrixXd::Identity(L,L));
    }

//    Eigen::SparseMatrix<Scalar> Id(L,L);
//    Id.setIdentity();
//    Id.makeCompressed();
//    MatrixType A_shift = (A_matrix - sigma * Eigen::MatrixXd::Identity(L,L));
    A_shift.makeCompressed();
    t_factorOp.toc();
    double t_factorOp_shift = t_factorOp.get_last_time_interval();
    t_factorOp.tic();
//    double sparcity_A_shift = (A_shift.toDense().array().cwiseAbs() > 1e-14 )
//                       .select(Eigen::MatrixXd::Ones(L,L),0).sum() / A_shift.size();
    lu_dense.compute(A_shift);
//    double sparcity_LU = (lu_dense.matrixLU().array().cwiseAbs() > 1e-14 )
//                                      .select(Eigen::MatrixXd::Ones(L,L),0).sum() / lu_dense.matrixLU().size();
//    std::cout << "sparcity_A_shift: " << sparcity_A_shift << std::endl;
//    std::cout << "sparcity_LU     : " << sparcity_LU      << std::endl;
//    std::cout << "LU              \n" << lu_dense.permutationP().inverse() * lu_dense.matrixLU()      << std::endl;
//    lu_sparse.compute(A_shift);
//    if (lu_sparse.info() == Eigen::ComputationInfo::Success){
//        sparseLU = true;
//    }else{
//        switch(lu_sparse.info()){
//            case Eigen::ComputationInfo::NoConvergence : {std::cerr << "NoConvergence: "    << lu_sparse.lastErrorMessage() << std::endl;}
//            case Eigen::ComputationInfo::InvalidInput  : {std::cerr << "InvalidInput: "     << lu_sparse.lastErrorMessage() << std::endl;}
//            case Eigen::ComputationInfo::NumericalIssue: {std::cerr << "NumericalIssue: "   << lu_sparse.lastErrorMessage() << std::endl;}
//            default :     {std::cout << "INFO: " << lu_sparse.info() << std::endl;}
//        }
//        lu_dense.compute(A_shift);
//        sparseLU = false;
//    }
    t_factorOp.toc();
    readyFactorOp = true;
    std::cout << "Time Factor Op [ms]: " << std::setprecision(3) << t_factorOp.get_last_time_interval() * 1000 << "  " << t_factorOp_shift * 1000 <<'\n';
//    std::cout << "Time Factor sh [ms]: " << std::setprecision(3)  '\n';
}

template<typename Scalar>
void SparseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
    Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
    switch (side){
        case Side::R: {
            if(sparseLU){x_out.noalias() = lu_sparse.solve(x_in);}
            else        {x_out.noalias() = lu_dense.solve(x_in);}
            break;
        }
        case Side::L: {
            std::cerr << "Left sided sparse shift invert hasn't been implemented yet..." << std::endl;
            exit(1);
            break;
        }
    }
    counter++;
}

template<typename Scalar>
void SparseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
    using namespace eigutils::eigSetting;
    switch (form){
        case Form::NONSYMMETRIC:
            switch (side) {
                case Side::R: {
                    Eigen::Map<VectorType> x_vec_in (x_in,  L);
                    Eigen::Map<VectorType> x_vec_out(x_out, L);
                    x_vec_out.noalias() = A_matrix * x_vec_in ;

                    break;
                }
                case Side::L: {
                    Eigen::Map<VectorTypeT> x_vec_in(x_in, L);
                    Eigen::Map<VectorTypeT> x_vec_out(x_out, L);
                    x_vec_out.noalias() = x_vec_in * A_matrix;
                    break;
                }
            }
            break;
        case Form::SYMMETRIC: {
            Eigen::Map<VectorType> x_vec_in (x_in,  L);
            Eigen::Map<VectorType> x_vec_out(x_out, L);
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
            break;
        }
    }
    counter++;
}



#ifdef EIGEN_USE_BLAS_SUSPEND
#define EIGEN_USE_BLAS
#undef EIGEN_USE_BLAS_SUSPEND
#endif




#endif //
