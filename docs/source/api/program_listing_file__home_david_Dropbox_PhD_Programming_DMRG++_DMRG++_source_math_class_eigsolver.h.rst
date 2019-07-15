
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_class_eigsolver.h:

Program Listing for File class_eigsolver.h
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_class_eigsolver.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/math/class_eigsolver.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-10-29.
   //
   
   #ifndef CLASS_EIGSOLVER_H
   #define CLASS_EIGSOLVER_H
   
   #include "math/arpack_extra/arpackpp_solver.h"
   #include "math/arpack_extra/matrix_recast.h"
   #include "nmspc_eigutils.h"
   
   template<typename T> class DenseMatrixProduct;
   template<typename T> class StlMatrixProduct;
   template<typename T> class SparseMatrixProduct;
   template<typename T> class DenseHamiltonianProduct;
   
   
   class class_eigsolver {
   private:
       size_t      logLevel  = 5;
   
   public:
   
       class_eigsolver(size_t logLevel_ = 5);
       void setLogLevel(size_t logLevelZeroToSix);
   
       eigutils::eigConfig     solverConf;
       eigutils::eigSolution   solution;
   
       int eig_dsyevd(const double* matrix, int L);
       int eig_dsyevd(double *matrix2eigvecs, double *eigvals, int L);
   
       int eig_zheevd(const std::complex<double>* matrix, int L);
       int eig_zheevd(std::complex<double>* matrix2eigvecs, double *eigvals, int L);
   
       int eig_dgeev (const double* matrix, int L);
       int eig_dgeev (const double* matrix, std::complex<double>* eigvecsR, std::complex<double>* eigvecsL, std::complex<double> *eigvals, int L);
       int eig_zgeev (const std::complex<double>* matrix, int L);
       int eig_zgeev (const std::complex<double>* matrix, std::complex<double>* eigvecsR, std::complex<double>* eigvecsL, std::complex<double> *eigvals, int L);
   
   
   
       template <typename Scalar>
       void subtract_phase(std::vector<Scalar> &eigvecs,int L, int nev)
       // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
       // is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
       // remove it from all other elements of v.
       {
           if constexpr (std::is_same<Scalar, std::complex<double>>::value) {
               if (nev > 0){
                   for (int i = 0; i < nev; i++) {
                       if (eigvecs[i*L].imag() == 0.0){continue;}
                       Scalar inv_phase = Scalar(0.0,-1.0) * std::arg(eigvecs[i * L]);
                       auto begin = eigvecs.begin() + i * L;
                       auto end = begin + L;
                       Scalar exp_inv_phase = std::exp(inv_phase);
                       std::transform(begin, end, begin,
                                      [exp_inv_phase](std::complex<double> num) -> std::complex<double>
                                      { return (num * exp_inv_phase); });
                       std::transform(begin, end, begin,
                                      [](std::complex<double> num) -> std::complex<double>
                                      { return std::abs(num.imag()) > 1e-15 ? num : std::real(num); });
                   }
               }else{
                   eigutils::eigLogger::log->error("Eigenvalues haven't been computed yet. Can't subtract phase.");
                   throw std::logic_error("Eigenvalues haven't been computed yet. Can't subtract phase.");
               }
           }
       }
   
      void eigs_init(const int L,
                      const int nev,
                      const int ncv,
                      const std::complex<double> sigma = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                      const eigutils::eigSetting::Type type = eigutils::eigSetting::Type::REAL,
                      const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
                      const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
                      const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                      const eigutils::eigSetting::Storage storage = eigutils::eigSetting::Storage::DENSE,
                      const bool compute_eigvecs_ = true,
                      const bool remove_phase_ = false
       );
   
       void eig_init(const eigutils::eigSetting::Type type        = eigutils::eigSetting::Type::REAL,
                     const eigutils::eigSetting::Form form        = eigutils::eigSetting::Form::SYMMETRIC,
                     const eigutils::eigSetting::Side side        = eigutils::eigSetting::Side::R,
                     const bool compute_eigvecs_                  = true,
                     const bool remove_phase_                     = false
       );
   
   
       template<eigutils::eigSetting::Type    type = eigutils::eigSetting::Type::REAL,
                eigutils::eigSetting::Form    form = eigutils::eigSetting::Form::SYMMETRIC,
                eigutils::eigSetting::Side    side = eigutils::eigSetting::Side::R,
                typename Derived>
       void eig(const Eigen::EigenBase<Derived> &matrix,
                const bool compute_eigvecs_           = true,
                const bool remove_phase_              = false);
   
       template<eigutils::eigSetting::Type    type = eigutils::eigSetting::Type::REAL,
                eigutils::eigSetting::Form    form = eigutils::eigSetting::Form::SYMMETRIC,
                eigutils::eigSetting::Side    side = eigutils::eigSetting::Side::R,
               typename Scalar>
       void eig(const Scalar * matrix,
                const int L,
                const bool compute_eigvecs_           = true,
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
       void eigs_dense(DenseHamiltonianProduct<Scalar> &matrix,
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
   
   
       template<typename Scalar>
       void eigs_stl(const Scalar *matrix,
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
       void eigs_stl(StlMatrixProduct<Scalar> &matrix,
                     const int nev,
                     const int ncv,
                     const std::complex<double> sigma      = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                     const eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC,
                     const eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::LM,
                     const eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
                     const bool compute_eigvecs_           = false,
                     const bool remove_phase_              = false,
                     Scalar *residual_                     = nullptr);
   
   
   
   };
   
   
   // Definitions
   
   
   
   template<eigutils::eigSetting::Type    type,
            eigutils::eigSetting::Form    form,
            eigutils::eigSetting::Side    side,
            typename Derived>
   void class_eigsolver::eig(const Eigen::EigenBase<Derived> &matrix,
                             const bool compute_eigvecs_,
                             const bool remove_phase_   )
   {
       eig<type,form>(matrix.derived().data(),matrix.rows(), compute_eigvecs_,remove_phase_);
   }
   
   
   
   
   
   template<eigutils::eigSetting::Type    type,
            eigutils::eigSetting::Form    form,
            eigutils::eigSetting::Side    side,
            typename Scalar>
   void class_eigsolver::eig(const Scalar * matrix,
                             const int L,
                             const bool compute_eigvecs_,
                             const bool remove_phase_   )
   {
       using namespace eigutils::eigSetting;
       eig_init(type,form,side,compute_eigvecs_,remove_phase_);
       int info = 0;
       try{
           if constexpr(form == Form::SYMMETRIC) {
               if constexpr(type == Type::REAL) {
                   static_assert(std::is_same<Scalar, double>::value);
                   info = eig_dsyevd(matrix,L);
               } else if constexpr (type == Type::CPLX) {
                   static_assert(std::is_same<Scalar, std::complex<double>>::value);
                   info = eig_zheevd(matrix,L);
               }
           }
   
           else
           if constexpr( form == Form::NONSYMMETRIC) {
               if constexpr(type == Type::REAL) {
                   static_assert(std::is_same<Scalar, double>::value);
                   info = eig_dgeev(matrix, L);
               } else if constexpr (type == Type::CPLX) {
                   static_assert(std::is_same<Scalar, std::complex<double>>::value);
                   info = eig_zgeev(matrix, L);
               }
           }
       }catch(std::exception &ex){
           eigutils::eigLogger::log->error("Eigenvalue solver failed: " + std::string(ex.what()) );
       }
   
       if (info == 0 and solverConf.remove_phase){
        // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
        // is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
        // remove it from all other elements of v.
           subtract_phase(solution.get_eigvecs<type,form,side>(),L, solution.meta.nev);
       }
   
   }
   
   
   
   
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
   
       eigs_init(L, nev, -1, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
   
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
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
       if constexpr(storage == Storage::DENSE){
           auto matrix_dense = DenseMatrixProduct<Scalar> (matrix,L);
           arpackpp_solver<DenseMatrixProduct<Scalar>> solver(matrix_dense, solverConf, solution,residual_);
           solver.eigs();
       }else if constexpr (storage == Storage::SPARSE){
           auto matrix_sparse = SparseMatrixProduct<Scalar> (matrix,L);
           arpackpp_solver<SparseMatrixProduct<Scalar>> solver(matrix_sparse, solverConf, solution,residual_);
           solver.eigs();
       }else if constexpr (storage == Storage::STL){
           auto matrix_stl = StlMatrixProduct<Scalar> (matrix,L);
           arpackpp_solver<StlMatrixProduct<Scalar>> solver(matrix_stl, solverConf, solution,residual_);
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
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
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
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
       arpackpp_solver<DenseMatrixProduct<Scalar>> solver(matrix_dense, solverConf, solution,residual_);
       solver.eigs();
   }
   
   template<typename Scalar>
   void class_eigsolver::eigs_dense      (DenseHamiltonianProduct<Scalar> &matrix,
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
       int L = matrix.rows();
       bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
       Type type = is_cplx ? Type::CPLX : Type::REAL;
       Storage storage = Storage::STL;
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
       arpackpp_solver<DenseHamiltonianProduct<Scalar>> solver(matrix, solverConf, solution,residual_);
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
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
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
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
       arpackpp_solver<SparseMatrixProduct<Scalar>> solver(matrix_sparse, solverConf, solution,residual_);
       solver.eigs();
   }
   
   
   template<typename Scalar>
   void class_eigsolver::eigs_stl   (const Scalar *matrix,
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
       Storage storage = Storage::STL;
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
       auto matrix_stl = SparseMatrixProduct<Scalar> (matrix,L);
       arpackpp_solver<StlMatrixProduct<Scalar>> solver(matrix_stl, solverConf, solution,residual_);
       solver.eigs();
   }
   
   
   
   
   template<typename Scalar>
   void class_eigsolver::eigs_stl       (StlMatrixProduct<Scalar> &matrix_stl,
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
       int L = matrix_stl.rows();
       bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
       Type type = is_cplx ? Type::CPLX : Type::REAL;
       Storage storage = Storage::STL;
       eigs_init(L, nev, ncv, sigma, type, form, ritz, side, storage, compute_eigvecs_, remove_phase_);
       arpackpp_solver<StlMatrixProduct<Scalar>> solver(matrix_stl, solverConf, solution,residual_);
       solver.eigs();
   }
   
   
   
   
   
   
   
   #endif //EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H
