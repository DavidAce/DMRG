//
// Created by david on 2018-05-06.
//

#ifndef CLASS_ARPACKPP_WRAPPER2_H
#define CLASS_ARPACKPP_WRAPPER2_H

#include <map>
#include <complex>
#include <vector>
#include <iostream>


namespace arpackpp{
    enum class Form{REAL_SYMMETRIC, REAL_GENERAL, COMPLEX_GENERAL};  // Real Symmetric, Real General or Complex General
//    enum class Side {L,R};           //Left or right eigenvectors
    enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE}; //choice of eigenvalue. LA is largest algebraic, and so on.
}
using namespace arpackpp;

//struct solverproperties{
//    Form form = Form::COMPLEX_GENERAL;
//    Side side = Side::R;
//    Ritz ritz = Ritz::LM;
//
//};


template<typename Scalar, Form form = Form::COMPLEX_GENERAL>
class class_arpackpp_wrapper2 {

private:
//    using Scalar = std::complex<double>;
    double  eigThreshold    = 1e-12;
    int     eigMaxIter      = 1000;
    int     ncv_max         = 10;
    int     Iter            = 0;
    int     rows            = 0;
    int     cols            = 0;
    int     iter            = 0;
    int     nev_found       = 0;
    int     n               = 0;
    int     counter         = 0;
    bool    eigvecs_found   = false;
    bool    eigvals_found   = false;
    std::vector<Scalar> eigvecs;
    std::vector<Scalar> eigvals;

    using  MapType = std::map<arpackpp::Ritz, std::string>;
    MapType RitzToString;
    template <typename Derived>
    void retrieve_solution(Derived & solution, int nev, bool compute_eigv) {
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
    public:



    class_arpackpp_wrapper2();

    const std::vector<Scalar> & get_eigvecs() const;
    const std::vector<Scalar> & get_eigvals() const;
    void dephase();

    void setThreshold(double newThreshold) {
        eigThreshold = newThreshold;
    }

    int GetIter(){
        return Iter;
    }

    void eig(Scalar *data,
             Ritz ritz,
             const int n,
             const int nev,
             const int ncv,
             bool compute_eigv=true);

    void eig(Scalar *data,
             std::string ritz,
             const int n,
             const int nev,
             const int ncv,
             bool compute_eigv=true);


};


#endif //DMRG_CLASS_ARPACKPP_WRAPPER2_H
