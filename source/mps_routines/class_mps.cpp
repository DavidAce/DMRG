//
// Created by david on 2017-11-13.
//
#include "class_mps.h"

void class_mps::initialize(const long local_dimension_){
    local_dimension = local_dimension_;
    GA.resize(array3{local_dimension,1,1});
    GA.setZero();
    GA(0,0,0) = 1;
    LA.resize(array1{1});
    LA.setConstant(1.0);

    GB.resize(array3{local_dimension,1,1});
    GB.setZero();
    GB(0,0,0) = 1;
    LB.resize(array1{1});
    LB.setConstant(1.0);
    swapped = false;
    L_tail = LB; /*! \todo , check?  was LA before */
}



Textra::Tensor<4,class_mps::Scalar> class_mps::get_theta() const {
    return asDiagonal(L_tail) //whatever L_A was in the previous step
            .contract(GA,             idx<1>({1},{1}))
            .contract(asDiagonal(LA), idx<1>({2},{0}))
            .contract(GB,             idx<1>({2},{1}))
            .contract(asDiagonal(LB), idx<1>({3},{0})).shuffle(array4{1,2,0,3});
    //Outputs:
    //      0  1
    //   2__|__|__3
}



void class_mps::swap_AB(){
    swapped    = !swapped;

    //Swap Gamma
    tmp3    = GA;
    GA      = GB;
    GB      = tmp3;

    //Swap Lambda
    tmp1    = LA;
    LA      = LB;
    LB      = tmp1;
    L_tail  = LB;

}