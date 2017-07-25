//
// Created by david on 7/21/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_MPS_H
#define FINITE_DMRG_EIGEN_CLASS_MPS_H

#include <n_tensor_extra.h>
#include <vector>
#include <iomanip>
#include <functional>
#include "funcs.h"

using namespace Textra;
using namespace std;



class Gamma_AB{
private:
    bool swap;
    Tensor3 tmp;
public:
    Tensor3 A,B;

    Gamma_AB(){
        swap = false;
    }
    void swap_AB(){
        swap = !swap;
        tmp = A;
        A = B;
        B = tmp;
    }
    //Access A and  B as 0 and 1 regardless of swap
    const Tensor3 &operator[](const long i)const{return i == 0 ? (swap ? B : A) : (swap? A : B);}
    Tensor3 &operator[](const long i){return i == 0 ?(swap ? B : A) : (swap? A : B);}


};

class Lambda_AB{
private:
    bool swap;
    Tensor1 tmp;
public:
    Tensor1 A, B;
    Lambda_AB(){
        swap = false;
    }
    void swap_AB() {
        swap = !swap;
        tmp = A;
        A = B;
        B = tmp;
    }
    //Access elements 0 and 1 regardless of swap
    const Tensor1 & operator[](const long i) const {return i == 0 ?(swap ? B : A) : (swap? A : B);}
    Tensor1 & operator[](const long i) {return i == 0 ?(swap ? B : A) : (swap? A : B);}
};



class class_MPS {
private:
    long local_dimension;
    long sites;
public:
    Gamma_AB  G;
    Lambda_AB L;
    Tensor1 L_tail;

    class_MPS(){
    }
    void initialize(const long dim, const long s){
        local_dimension = dim;
        sites = s;
        for(unsigned int site = 0; site < sites; site++){
            G[site].resize(array3{local_dimension,1,1});
            G[site].setZero();
            G[site](0,0,0) = 1;
            L[site].resize(array1{1});
            L[site].setConstant(1.0);
            L_tail = L[site];
        }
    }

    Tensor4 get_theta() const {
        return asDiagonal(L_tail) //whatever L_A was in the previous step
                .contract(G.A, idxlist1{idx2(1,1)})
                .contract(asDiagonal(L.A), idxlist1{idx2(2,0)})
                .contract(G.B, idxlist1{idx2(2,1)})
                .contract(asDiagonal(L.B), idxlist1{idx2(3,0)})
                .shuffle(array4{1,0,2,3});
    }

    void swap_AB(){
        L_tail = L.A;
        G.swap_AB();
        L.swap_AB();
    }
    void print_error(const Tensor4 &Hamiltonian) {
        std::cout << std::setprecision(16);
//        Eigen::Array<Scalar,2, 1> E;
//        for(long site = 0; site < sites; site++) {
            Tensor4 LGLGL = asDiagonal(L_tail)
                    .contract(G.A, idxlist1{idx2(1,1)})
                    .contract(asDiagonal(L.A), idxlist1{idx2(2,0)})
                    .contract(G.B, idxlist1{idx2(2,1)})
                    .contract(asDiagonal(L.B), idxlist1{idx2(3,0)})
                    .shuffle(array4{1,2,0,3});

            Tensor0 result = Hamiltonian.contract(LGLGL.conjugate(), idxlist2{idx2(0,0), idx2(1,1)})
                    .contract(LGLGL, idxlist4{idx2(0,0),idx2(1,1),idx2(2,2),idx2(3,3)});
//            swap_AB();
//            E(site) = result(0);
//        }
        cout << "E_iDMRG = " << result(0) << endl;
        cout << "E_exact = " << -1.063544409973372 << endl;
        cout << std::scientific;
        cout << "Error   = " << -1.063544409973372 - result(0) << endl;

    }



};


#endif //FINITE_DMRG_EIGEN_CLASS_MPS_H
