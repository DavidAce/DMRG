//
// Created by david on 2019-07-06.
//

#ifndef DMRG_CLASS_VIDAL_SITE_H
#define DMRG_CLASS_VIDAL_SITE_H

#include <general/nmspc_tensor_extra.h>



class class_vidal_site {
public:
    using Scalar = std::complex<double>;

private:
    std::optional<size_t> position;
    Eigen::Tensor<Scalar,3> G;                  /*!< \f$\Gamma \f$*/
    Eigen::Tensor<Scalar,1> L;                  /*!< \f$\Lambda\f$*/
public:

    class_vidal_site() = default;
    class_vidal_site(
            const Eigen::Tensor<Scalar,3> &G_,
            const Eigen::Tensor<Scalar,1> &L_,
            size_t pos):
            position(pos),G(G_),L(L_){};


    bool isReal()const;
    const Eigen::Tensor<Scalar,3> &get_G() const;
    const Eigen::Tensor<Scalar,1> &get_L() const;
    Eigen::Tensor<Scalar,3> &get_G();
    Eigen::Tensor<Scalar,1> &get_L();

    Eigen::Tensor<Scalar,3> get_A()  const;
    Eigen::Tensor<Scalar,3> get_B()  const;

    std::tuple<long,long,long> get_dims() const;

    long get_spin_dim() const;
    long get_chiL()     const;
    long get_chiR()     const;

    void set_position(const size_t position_);
    size_t get_position() const;

    void set_mps(const Eigen::Tensor<Scalar,3> &G_, const Eigen::Tensor<Scalar,1> &L_);
    void set_L(const Eigen::Tensor<Scalar,1> &L_);
    void set_G(const Eigen::Tensor<Scalar,3> &G_);





};


#endif //DMRG_CLASS_VIDAL_SITE_H
