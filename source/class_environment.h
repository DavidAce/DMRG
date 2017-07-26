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
    std::string picture;
    int size;         //Size is the number of particles that have been contracted into this environment

    Tensor3 block;
    class_environment_L();
    void enlarge(const class_MPS &MPS, const Tensor4 &W);
};


class class_environment_R {
private:

public:
    std::string single_picture = "-";
    std::string picture;
    int size;         //Size is the number of particles that have been contracted into this environment

    Tensor3 block;
    class_environment_R();

    void enlarge(const class_MPS &MPS,const Tensor4 &W);
};


#endif //FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H
