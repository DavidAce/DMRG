//
// Created by david on 2018-05-06.
//

#ifndef CLASS_ARPACK_EIGSOLVER_H
#define CLASS_ARPACK_EIGSOLVER_H

#include <map>
#include <complex>
#include <vector>
#include <iostream>
#include <memory>
#include <general/class_tic_toc.h>
#include "nmspc_eigsolver_props.h"


class class_superblock;

template<typename Scalar, Form form = Form::GENERAL>
class class_eigsolver_arpack {

private:
//    using T = std::complex<double>;
    double  eigThreshold    = 1e-12;
    int     eigMaxIter      = 1000;
    int     eigMaxNcv       = 16;
    int     Iter            = 0;
    int     rows            = 0;
    int     cols            = 0;
    int     iter            = 0;
    int     nev_found       = 0; // Found eigenvectors. aka cols.
    int     n               = 0; // Linear dimension of the input matrix to diagonalize, aka rows.
    int     counter         = 0;
    bool    eigvecs_found   = false;
    bool    eigvals_found   = false;
    bool    compute_eigvecs = false;
    bool    remove_phase    = false;
    Scalar  *residual = nullptr;
    std::vector<Scalar> eigvecs;
    std::vector<Scalar> eigvals;

    using  MapType = std::map<eigsolver_properties::Ritz, std::string>;
    MapType RitzToString;
    char ritz_char[2];

    template <typename Derived>
    void find_solution(Derived &solution, int nev, bool compute_eigv) {
        if (compute_eigv) {
            solution.FindEigenvectors();
            eigvals_found = solution.EigenvaluesFound();
            eigvecs_found = solution.EigenvectorsFound();
            Iter = solution.GetIter();
            n = solution.GetN();
            nev_found = std::min(nev, solution.GetNev());
            rows = n;
            cols = nev_found;
            eigvecs.insert(eigvecs.begin(), solution.RawEigenvectors(), solution.RawEigenvectors() + n * nev_found);
            eigvals.insert(eigvals.begin(), solution.RawEigenvalues(), solution.RawEigenvalues() + nev_found);
        }else{
            solution.FindEigenvalues();
            eigvals_found = solution.EigenvaluesFound();
            Iter = solution.GetIter();
            n = solution.GetN();
            nev_found = std::min(nev, solution.GetNev());
            rows = n;
            cols = nev_found;
            int cols = std::min(nev, solution.GetNev());
            eigvals.insert(eigvals.begin(), solution.RawEigenvalues(), solution.RawEigenvalues() + cols);
        }
    }
    void shift_invert_eigvals(Scalar sigma);
    void subtract_phase();

    class_tic_toc t_sol;
    class_tic_toc t_get;
    class_tic_toc t_sub;
    class_tic_toc t_all;

public:



    class_eigsolver_arpack();


    class_eigsolver_arpack(
            const double eigThreshold_,
            const int eigMaxIter_,
            const int eigMaxNcv_,
            const bool get_eigvecs_=false,
            const bool remove_phase_=false);


    const std::vector<Scalar> & ref_eigvecs() const;
    const std::vector<Scalar> & ref_eigvals() const;
    auto ref_eig_vecs_vals() const{
        return std::make_pair(ref_eigvecs(), ref_eigvals());
    }

    const std::vector<Scalar> get_eigvecs() const;
    const std::vector<Scalar> get_eigvals() const;
    auto get_eig_vecs_vals() const{
        return std::make_pair(get_eigvecs(), get_eigvals());
    }



    void setThreshold(double newThreshold) {
        eigThreshold = newThreshold;
    }
    void setMaxIter(int newMaxIter) {
        eigMaxIter = newMaxIter;
    }

    void setMaxNcv(int newMaxNcv) {
        eigMaxNcv = newMaxNcv;
    }


    int GetIter(){
        return Iter;
    }

    int GetNevFound(){
        return nev_found;
    }

    int Rows(){
        return rows;
    }

    int Cols(){
        return cols;
    }

    void eig(const Scalar *matrix_data,
             const int n,
             const int nev,
             const int ncv,
             const Ritz ritz = Ritz::SR,
             const Side side = Side::R,
             const bool compute_eigvecs_=false,
             const bool remove_phase_=false,
             Scalar *residual_ = nullptr
    );


    void eig_shift_invert(
             Scalar *matrix_data,
             const int n,
             const int nev,
             const int ncv,
             const Scalar shift,
             const Ritz ritz,
             const bool compute_eigvecs_= false,
             const bool remove_phase_= false,
             Scalar *residual_ = nullptr
    );



    void eig_shift_invert2(
            const Scalar *Lblock_,                   /*!< The left block tensor.  */
            const Scalar *Rblock_,                   /*!< The right block tensor.  */
            const Scalar *HA_,                       /*!< The left Hamiltonian MPO's  */
            const Scalar *HB_,                       /*!< The right Hamiltonian MPO's */
            const std::array<long,4> shape_theta4_,  /*!< An array containing the shapes of theta  */
            const std::array<long,4> shape_mpo4_ ,   /*!< An array containing the shapes of the MPO  */
            const int nev,
            const int ncv,
            const Scalar shift,
            const Ritz ritz,
            const bool compute_eigvecs_= false,
            const bool remove_phase_= false,
            Scalar *residual_ = nullptr
    );


    void optimize_mps(
            const Scalar *Lblock_,                   /*!< The left block tensor.  */
            const Scalar *Rblock_,                   /*!< The right block tensor.  */
            const Scalar *HA_,                       /*!< The left Hamiltonian MPO's  */
            const Scalar *HB_,                       /*!< The right Hamiltonian MPO's */
            const std::array<long,4> shape_theta4_,  /*!< An array containing the shapes of theta  */
            const std::array<long,4> shape_mpo4_ ,   /*!< An array containing the shapes of the MPO  */
            const int nev,
            const int ncv,
            const Ritz ritz = Ritz::SR,
            const bool remove_phase_ = true,
            Scalar *residual_ = nullptr);

};


#endif //CLASS_ARPACK_EIGSOLVER_H