//
// Created by david on 2018-10-30.
//

#ifndef CLASS_EIGSOLVER_BASE_H
#define CLASS_EIGSOLVER_BASE_H
#include <complex>
#include <vector>
#include <map>
#include <memory>
#include <general/class_tic_toc.h>
#include "general/nmspc_eigutils.h"
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


    void shift_invert_eigvals(std::complex<double> sigma);
    void subtract_phase();
    void eigs();
    template <typename Derived>  void find_solution(Derived &solver, int nev);
    template <typename Derived>  void copy_solution_symm(Derived &solver);
    template <typename Derived>  void copy_solution_nsym(Derived &solver);


    MatrixType               &matrix;
    eigutils::eigConfig      &solverConf;
    eigutils::eigSolution    &solution;

    arpackpp_solver(
            MatrixType              &matrix_,
            eigutils::eigConfig     &solverConf_,
            eigutils::eigSolution   &solution_);
};









//Definitions












#endif //EIGBENCH_CLASS_EIGSOLVER_BASE_H
