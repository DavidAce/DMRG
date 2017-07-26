#ifndef FINITE_DMRG_EIGEN_CLASS_MPS_H
#define FINITE_DMRG_EIGEN_CLASS_MPS_H

#include <n_tensor_extra.h>

using namespace Textra;
using namespace std;



class Gamma_AB{
private:
    bool swap;
    Tensor3 tmp;
public:
    Tensor3 A,B;
    Gamma_AB(){swap = false;}
    void swap_AB(){
        swap = !swap;
        tmp = A;
        A = B;
        B = tmp;
    }
    //Access A and  B as 0 and 1 regardless of swap
    const Tensor3 &operator[](const int i)const{return i == 0 ? (swap ? B : A) : (swap? A : B);}
    Tensor3 &operator[](const int i){return i == 0 ?(swap ? B : A) : (swap? A : B);}


};

class Lambda_AB{
private:
    bool swap;
    Tensor1 tmp;
public:
    Tensor1 A, B;
    Lambda_AB(){swap = false;}
    void swap_AB() {
        swap = !swap;
        tmp = A;
        A = B;
        B = tmp;
    }
    //Access elements 0 and 1 regardless of swap
    const Tensor1 & operator[](const int i) const {return i == 0 ?(swap ? B : A) : (swap? A : B);}
    Tensor1 & operator[](const int i) {return i == 0 ?(swap ? B : A) : (swap? A : B);}
};



class class_MPS {
private:
    long local_dimension;
    long sites;
public:
    Gamma_AB  G;
    Lambda_AB L;
    Tensor1 L_tail;

    class_MPS(){}
    void initialize(const long dim, const long s);
    Tensor4 get_theta() const;
    void swap_AB();
    void print_error(const Tensor4 &Hamiltonian);



};


#endif //FINITE_DMRG_EIGEN_CLASS_MPS_H
