
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_arpack_extra_arpackpp_solver.h:

Program Listing for File arpackpp_solver.h
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_arpack_extra_arpackpp_solver.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/math/arpack_extra/arpackpp_solver.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-10-30.
   //
   
   #ifndef CLASS_ARPACKPP_SOLVER_H
   #define CLASS_ARPACKPP_SOLVER_H
   #include <vector>
   #include <general/class_tic_toc.h>
   #include "math/nmspc_eigutils.h"
   #include <complex>
   
   #define profile_arpack 0
   
   
   template<typename MatrixType>
   class arpackpp_solver {
   private:
       int nev_internal;
       int ncv_internal;
   
   
       void eigs_sym();
       void eigs_nsym();
       void eigs_comp();
   
   public:
   
       using Scalar = typename MatrixType::Scalar;
       class_tic_toc t_sol;
       class_tic_toc t_get;
       class_tic_toc t_sub;
       class_tic_toc t_all;
   
   
   //    void shift_invert_eigvals(std::complex<double> sigma);
   //    void subtract_phase(std::vector<Scalar> & eigvecs);
       void eigs();
       template <typename Derived>  void find_solution(Derived &solver, int nev);
   
   //    template <typename Derived>  void copy_solution_symm(Derived &solver);
   //    template <typename Derived>  void copy_solution_nsym(Derived &solver);
   
   //    template <eigutils::eigSetting::Type type,
   //              eigutils::eigSetting::Form form,
   //              eigutils::eigSetting::Side side,
   //              typename Derived>
   //    void copy_solution(Derived &solver);
   
   
       MatrixType               &matrix;
       eigutils::eigConfig      &solverConf;
       eigutils::eigSolution    &solution;
       Scalar                   *residual;
       arpackpp_solver(
               MatrixType              &matrix_,
               eigutils::eigConfig     &solverConf_,
               eigutils::eigSolution   &solution_,
               Scalar                  *residual_
               );
   
   
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
   
   
       template <eigutils::eigSetting::Type type,
               eigutils::eigSetting::Form form,
               eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R,
               typename Derived>
       void copy_solution(Derived &solver){
           using namespace eigutils::eigSetting;
           auto & eigvecs = solution.get_eigvecs<type,form,side>();
           auto & eigvals = solution.get_eigvals<form>();
           int eigvecsize = solution.meta.rows * solution.meta.cols;
           int eigvalsize = solution.meta.cols;
           solution.meta.form = form;
           solution.meta.type = type;
           solution.meta.side = side;
           // Copy eigenvalues
           eigvals.resize(eigvalsize);
           if constexpr(form == Form::SYMMETRIC){
               if constexpr (type == Type::REAL){
                   std::copy(solver.RawEigenvalues() ,solver.RawEigenvalues() + eigvalsize,eigvals.begin());
               }else {
                   for (int j = 0; j < solution.meta.cols; j++) {
   //                    assert(std::abs(solver.Eigenvalue(j).imag()) < 1e-15 and "Discarding imaginary part!" );
                       if(std::abs(solver.Eigenvalue(j).imag()) > solverConf.eigThreshold){
                           std::cerr << "WARNING: Discarding imaginary part: " << solver.Eigenvalue(j).imag() << std::endl;
                       }
                       eigvals[j] = solver.Eigenvalue(j).real();
                   }
               }
           }else if constexpr(form == Form::NONSYMMETRIC){
               if constexpr(type == Type::REAL){
                   for (int j = 0; j < solution.meta.cols; j++) {
                       eigvals[j] = std::complex<double>(solver.EigenvalueReal(j), solver.EigenvalueImag(j));
                   }
               }
               else if constexpr(type == Type::CPLX){
                   for (int j = 0; j < solution.meta.cols; j++) {
                       eigvals[j] = solver.Eigenvalue(j);
                   }
               }
           }
   
   
   
           // Copy eigenvectors
           if (solverConf.compute_eigvecs){
               eigvecs.resize(eigvecsize);
               if constexpr (type == Type::REAL and form == Form::SYMMETRIC){
                       std::copy(solver.RawEigenvectors(),solver.RawEigenvectors() + eigvecsize, eigvecs.begin());
               }else{
                   int count = 0;
                   for (int j = 0; j < solution.meta.cols; j++){
                       for (int i = 0; i < solution.meta.rows; i++){
                           eigvecs[count++] = solver.Eigenvector(j,i);
                       }
                   }
               }
          }
   
           if(solverConf.remove_phase){
               subtract_phase(eigvecs,solution.meta.rows,solution.meta.cols);
           }
       }
   
   
   
   };
   
   
   
   
   
   
   
   
   
   //Definitions
   
   
   
   
   
   
   
   
   
   
   
   
   #endif //CLASS_ARPACKPP_SOLVER_H
