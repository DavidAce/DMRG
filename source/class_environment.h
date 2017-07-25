//
// Created by david on 7/21/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H
#define FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H

#include <n_tensor_extra.h>
#include "class_MPS.h"


//enum class environment_side{LEFT,RIGHT};

using namespace Textra;
using namespace std;
class class_environment_L {
private:

public:
    std::string single_picture = "=";
    std::string full_picture;
    size_t size;         //Size is the number of particles that have been contracted into this environment

    Tensor3 block;
    class_environment_L();
//    class_environment_L& operator=(class_environment_L other){
//        std::swap(*this, other); // (2)
//
//        return *this;
//    }
    // move constructor
//    class_environment_L(class_environment_L&& other): class_environment_L() // initialize via default constructor, C++11 only
//    {
//        std::swap(*this, other);
//    }



    void enlarge(const class_MPS &MPS, const Tensor4 &W);

//    void enlarge_left(const Gamma_AB &G, const Lambda_AB &L, const Tensor4 &W);
//    void enlarge_right(const Gamma_AB &G, const Lambda_AB &L, const Tensor4 &W);

};


class class_environment_R {
private:

public:
    std::string single_picture = "-";
    std::string full_picture;
    size_t size;         //Size is the number of particles that have been contracted into this environment

    Tensor3 block;
    class_environment_R();


//    class_environment_R& operator=(class_environment_R other){
//        std::swap(*this, other); // (2)
//
//        return *this;
//    }
    // move constructor
//    class_environment_R(class_environment_R&& other): class_environment_R() // initialize via default constructor, C++11 only
//    {
//        std::swap(*this, other);
//    }


    void enlarge(const class_MPS &MPS,const Tensor4 &W);

//    void enlarge_left(const Gamma_AB &G, const Lambda_AB &L, const Tensor4 &W);
//    void enlarge_right(const Gamma_AB &G, const Lambda_AB &L, const Tensor4 &W);

};


#endif //FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H
