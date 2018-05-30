//
// Created by david on 2018-05-06.
//

#ifndef CLASS_ARPACKPP_WRAPPER2_H
#define CLASS_ARPACKPP_WRAPPER2_H

#include <map>
#include <complex>
#include <vector>
#include <iostream>
#include "class_tic_toc.h"

namespace arpackpp{
    enum class Form{SYMMETRIC, GENERAL};  // Real Symmetric, Real General or Complex General
    enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE}; //choice of eigenvalue. C is largest algebraic, and so on.
    enum class Side {L,R};           //Left or right eigenvectors
}

using namespace arpackpp;

template<typename Scalar, Form form = Form::GENERAL>
class class_arpack_eigsolver {

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
    std::vector<Scalar> eigvecs;
    std::vector<Scalar> eigvals;

    using  MapType = std::map<arpackpp::Ritz, std::string>;
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

    class_tic_toc t_sol;
    class_tic_toc t_get;
    class_tic_toc t_sub;
    class_tic_toc t_all;

public:



    class_arpack_eigsolver();
    class_arpack_eigsolver(const Scalar *matrix_data,
                           const Ritz ritz,
                           const Side side,
                           const int n,
                           const int nev,
                           const int ncv,
                           const double threshold,
                           const int MaxIter,
                           const bool getvecs=false,
                           const bool dephase=false,
                           Scalar *residp = NULL);
    class_arpack_eigsolver(const Scalar *matrix_data,
                           const Ritz ritz,
                           const Side side,
                           const int n,
                           const int nev,
                           const int ncv,
                           const bool getvecs=false,
                           const bool dephase=false,
                           Scalar *residp = NULL);


    class_arpack_eigsolver(const Scalar *Lblock,                            /*!< The left block tensor.  */
                           const Scalar *Rblock,                            /*!< The right block tensor.  */
                           const Scalar *HA,                                /*!< The left Hamiltonian MPO's  */
                           const Scalar *HB,                                /*!< The right Hamiltonian MPO's */
                           const std::array<long,4> shape_theta4,           /*!< An array containing the shapes of theta  */
                           const std::array<long,4> shape_mpo4 ,            /*!< An array containing the shapes of the MPO  */
                           const Ritz ritz,
                           const int nev,
                           const int ncv,
                           const bool bool_dephase=true,
                           Scalar *resid = nullptr);

    class_arpack_eigsolver(const Scalar *Lblock,        /*!< The left block tensor.  */
                           const Scalar *Rblock,        /*!< The right block tensor.  */
                           const Scalar *HA,            /*!< The left Hamiltonian MPO's  */
                           const Scalar *HB,            /*!< The right Hamiltonian MPO's */
                           const std::array<long,4> shape_theta4,         /*!< An array containing the shapes of theta  */
                           const std::array<long,4> shape_mpo4 ,           /*!< An array containing the shapes of the MPO  */
                           const Ritz ritz,
                           const int nev,
                           const int ncv,
                           const double threshold,
                           const int MaxIter,
                           const bool bool_dephase=true,
                           Scalar *resid = nullptr);




    const std::vector<Scalar> & ref_eigvecs() const;
    const std::vector<Scalar> & ref_eigvals() const;

    const std::vector<Scalar> get_eigvecs() const;
    const std::vector<Scalar> get_eigvals() const;


//

    void subtract_phase();

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

    void eig(const Scalar *matrix_data,
             const Ritz ritz,
             const Side side,
             const int n,
             const int nev,
             const int ncv,
             const double threshold,
             const int  MaxIter,
             const bool bool_find_eigvecs=false,
             const bool bool_dephase=false,
             Scalar *residp = nullptr
    );

    const std::pair<const std::vector<Scalar>&,const std::vector<Scalar>&>
    eig_ref_vec_val(const Scalar *matrix_data,
                    const Ritz ritz,
                    const Side side,
                    const int n,
                    const int nev,
                    const int ncv,
                    const bool bool_dephase = false,
                    Scalar *residp = nullptr
    );

    const std::pair<const std::vector<Scalar>,const std::vector<Scalar>>
    eig_get_vec_val(const Scalar *matrix_data,
                    const Ritz ritz,
                    const Side side,
                    const int n,
                    const int nev,
                    const int ncv,
                    const bool bool_dephase = false,
                    Scalar *residp = nullptr
    );

    void optimize_mps(
            const Scalar *Lblock_,        /*!< The left block tensor.  */
            const Scalar *Rblock_,        /*!< The right block tensor.  */
            const Scalar *HA_,            /*!< The left Hamiltonian MPO's  */
            const Scalar *HB_,            /*!< The right Hamiltonian MPO's */
            const std::array<long,4> shape_theta4_,         /*!< An array containing the shapes of theta  */
            const std::array<long,4> shape_mpo4_ ,           /*!< An array containing the shapes of the MPO  */
            const Ritz ritz,
            const int nev,
            const int ncv,
            const double threshold,
            const int MaxIter,
            const bool bool_dephase = true,
            Scalar *residp = nullptr);

};


#endif //DMRG_CLASS_ARPACKPP_WRAPPER2_H