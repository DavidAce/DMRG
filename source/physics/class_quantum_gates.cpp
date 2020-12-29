//
// Created by david on 2020-12-20.
//

#include "class_quantum_gates.h"
#include "general/nmspc_tensor_extra.h"
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

qm::Gate::Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_):op(op_),pos(pos_){}
qm::Gate::Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_, Scalar delta):pos(pos_){
    op = Textra::MatrixToTensor((delta * Textra::TensorMatrixMap(op_)).exp());
}

void qm::Gate::exp_inplace(Scalar delta) { op = Textra::MatrixToTensor((delta * Textra::TensorMatrixMap(op)).exp());}
qm::Gate qm::Gate::exp(Scalar delta) const { return Gate(op,pos,delta);}
qm::Gate qm::Gate::adjoint() const { return Gate(op.conjugate().shuffle(Textra::array2{1, 0}), pos); }
bool qm::Gate::isUnitary(double prec) const { return Textra::TensorMatrixMap(op).isUnitary(prec); }
