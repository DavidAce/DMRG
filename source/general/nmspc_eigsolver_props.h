//
// Created by david on 2018-06-07.
//

#ifndef NMSPC_EIGSOLVER_PROPS_H
#define NMSPC_EIGSOLVER_PROPS_H

namespace eigsolver_properties{
    enum class Form{SYMMETRIC, GENERAL};  // Real Symmetric, Real General or Complex General
    enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE}; //choice of eigenvalue. C is largest algebraic, and so on.
    enum class Side {L,R};           //Left or right eigenvectors
}

//using namespace eigsolver_properties;
#endif //NMSPC_EIGSOLVER_PROPS_H
