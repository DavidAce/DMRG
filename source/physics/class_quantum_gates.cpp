//
// Created by david on 2020-12-20.
//

#include "class_quantum_gates.h"
#include "general/nmspc_tensor_extra.h"
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

Eigen::Tensor<qm::Scalar,2> qm::Gate::exp_internal(const Eigen::Tensor<Scalar,2> & op_, Scalar alpha) const {
    /* Note fore flbit:
    *  Let h = op(i,i), i.e. h are the diagonal entries in op
    *  Let alpha = -i * delta be purely imaginary (since delta is real). So delta = imag(-alpha)
    *  Now notice that
    *       exp( -i * delta * h )
    *  can become imprecise when delta is large, since delta and h are real and h ~O(1)
    *
    *  Instead, in the flbit case, we can exploit that h is real to compute
    *       exp(-i * mod(delta * real(h), 2*pi))
    *  which is equivalent but with a much smaller exponent.
    *
    *  Remember to do the modulo separately on each diagonal entry h!
    */
    auto op_map = Textra::TensorMatrixMap(op_);
    if(op_map.isDiagonal() and op_map.imag().isZero() and std::real(alpha) == 0){
        auto minus_i = std::complex<double>(0,-1);

        auto diag = op_map
            .diagonal()
            .unaryViewExpr([&alpha, &minus_i](const Scalar & h)
                           {return std::exp(minus_i * std::fmod(std::imag(-alpha)*std::real(h),2.0*M_PI));})
            .asDiagonal();
        return Textra::MatrixToTensor(diag);
//        std::cout << "alpha: " << alpha << std::endl;
//        std::cout << "op:\n" << op << std::endl;
//        std::cout << "old:\n" << (alpha * Textra::TensorMatrixMap(op_)).exp() << std::endl;
//        if(not Textra::TensorMatrixMap(op).isApprox((alpha * Textra::TensorMatrixMap(op_)).exp())){
//            throw std::runtime_error("Mismatch");
//        }
    }else{
        return Textra::MatrixToTensor((alpha * Textra::TensorMatrixMap(op_)).exp());
    }


}


qm::Gate::Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_):op(op_),pos(pos_){}
qm::Gate::Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_, Scalar alpha):pos(pos_){
    op = exp_internal(op_,alpha);
}

void qm::Gate::exp_inplace(Scalar alpha) { op = exp_internal(op,alpha);}
qm::Gate qm::Gate::exp(Scalar alpha) const { return Gate(op,pos,alpha);}
bool qm::Gate::isUnitary(double prec) const { return Textra::TensorMatrixMap(op).isUnitary(prec); }
Eigen::Tensor<qm::Scalar,2> & qm::Gate::adjoint() const {
    if(adj) return adj.value();
    adj = op.conjugate().shuffle(Textra::array2{1, 0});
    return adj.value();
}
