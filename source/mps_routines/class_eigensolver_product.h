//
// Created by david on 2017-11-15.
//

#ifndef DMRG_CLASS_GROUNDSTATE_SEARCH_H
#define DMRG_CLASS_GROUNDSTATE_SEARCH_H


#include "class_superblock.h"

class class_eigensolver_product {
private:
    using Scalar = double;
    const Textra::Tensor<3,Scalar> &Lblock;
    const Textra::Tensor<3,Scalar> &Rblock;
    const Textra::Tensor<6,Scalar> &WW;
    const Textra::array<4>         &shape4;
    Textra::array<1> shape1;
public:
    explicit class_eigensolver_product(const Textra::Tensor<3,Scalar> &Lblock_,
                                       const Textra::Tensor<3,Scalar> &Rblock_,
                                       const Textra::Tensor<6,Scalar> &WW_,
                                       const Textra::array<4> &shape4_):
            Lblock(Lblock_),
            Rblock(Rblock_),
            WW(WW_),
            shape4(shape4_)
    {
        shape1[0] = shape4[0]*shape4[1]*shape4[2]*shape4[3];
    };

//Functions for eigenvalue solver Spectra
    int rows()const {return (int)shape1[0];}
    int cols()const {return (int)shape1[0];}
    void perform_op(const double *x_in, double *y_out) const {
        Eigen::TensorMap<Textra::const_Tensor<4,double>> x_input (x_in, shape4);
        Eigen::TensorMap<Textra::Tensor<1,double>> y_output (y_out, shape1);
        y_output = Lblock
                        .contract(x_input,   idx<1>({1},{1}))
                        .contract(WW ,       idx<3>({1,2,3},{0,4,5}))
                        .contract(Rblock,    idx<2>({1,2},{1,2}))
                        .shuffle(array4{1,0,2,3})
                        .reshape(shape1);
//        Textra::Tensor<1,double> y_output =
//                Lblock
//                .contract(x_input,   idx<1>({1},{1}))
//                .contract(WW ,       idx<3>({1,2,3},{0,4,5}))
//                .contract(Rblock,    idx<2>({1,2},{1,2}))
//                .shuffle(array4{1,0,2,3})
//                .reshape(shape1);
//        std::move(y_output.data(), y_output.data()+y_output.size(),  y_out);
    }
};


#endif //DMRG_CLASS_GROUNDSTATE_SEARCH_H
