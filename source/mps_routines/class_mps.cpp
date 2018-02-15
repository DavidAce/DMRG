//
// Created by david on 2017-11-13.
//
#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif
#include "class_mps.h"

using namespace std;
using namespace Textra;



// Function definitions

template<typename Scalar>
void class_mps<Scalar>::initialize(const long local_dimension_){
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

template<typename Scalar>
void class_mps<Scalar>::swap_AB(){
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

template<typename Scalar>
Tensor<Scalar,3> class_mps<Scalar>::A() const{
    return asDiagonal(L_tail).contract(GA, idx<1>({1},{1})).shuffle(array3{1,0,2});
};
template<typename Scalar>
Tensor<Scalar,3> class_mps<Scalar>::B() const{
    return GB.contract(asDiagonal(LB), idx<1>({2},{0}));
};

template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::thetaL() const{
    return A()
            .contract(asDiagonal(LA),idx<1>({2},{0}) )
            .contract(GB, idx<1>({2},{1}))
            .shuffle(array4{0,2,1,3});
};

template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::thetaR() const{
    return GA
            .contract(asDiagonal(LA),idx<1>({2},{0}) )
            .contract(B(), idx<1>({2},{1}))
            .shuffle(array4{0,2,1,3});
};


template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::get_theta() const {
    return asDiagonal(L_tail) //whatever L_A was in the previous step
            .contract(GA,             idx<1>({1},{1}))
            .contract(asDiagonal(LA), idx<1>({2},{0}))
            .contract(GB,             idx<1>({2},{1}))
            .contract(asDiagonal(LB), idx<1>({3},{0})).shuffle(array4{1,2,0,3});

    //Outputs:
    //      0  1
    //   2__|__|__3
}

template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::get_transfer_matrix_L()const{
    return A().contract(A().conjugate(), idx<1>({0},{0}))
              .shuffle(array4{0,2,1,3});
};

template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::get_transfer_matrix_R()const{
    return B().contract(B().conjugate(), idx<1>({0},{0}))
              .shuffle(array4{0,2,1,3});;
};


template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::get_transfer_2_site_matrix_L()const{
    return thetaL().contract(thetaL().conjugate(), idx<2>({0,1},{0,1})).shuffle(array4{0,2,1,3});
};


template<typename Scalar>
Tensor<Scalar,4> class_mps<Scalar>::get_transfer_2_site_matrix_R()const{
    return thetaR().contract(thetaR().conjugate(), idx<2>({0,1},{0,1})).shuffle(array4{0,2,1,3});
};

// Explicit instantiations

template class class_mps<double>;
template class class_mps<std::complex<double>>;

