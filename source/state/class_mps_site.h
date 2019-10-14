//
// Created by david on 2019-07-06.
//

#ifndef DMRG_CLASS_MPS_SITE_H
#define DMRG_CLASS_MPS_SITE_H

#include <general/nmspc_tensor_extra.h>



class class_mps_site {
public:
    using Scalar = std::complex<double>;
private:
    Eigen::Tensor<Scalar,3> M;                  /*!< \f$M\f$ A or B tensor*/
    Eigen::Tensor<Scalar,1> L;                  /*!< \f$\Lambda\f$*/
    std::optional<size_t> position;
    std::optional<Eigen::Tensor<Scalar,1>> LC;  /*!< \f$\Lambda_C\f$ Center lambda, if this is a center matrix*/
    mutable std::optional<Eigen::Tensor<Scalar,3>> MC;
    //    Eigen::Tensor<Scalar,3> G;                  /*!< \f$\Gamma \f$*/

public:

    class_mps_site() = default;
    class_mps_site(
            const Eigen::Tensor<Scalar,3> &M_,
            const Eigen::Tensor<Scalar,1> &L_,
            size_t pos)
            :
            M(M_),
            L(L_),
            position(pos){};



    [[nodiscard]] bool isReal() const;
    [[nodiscard]] bool isCenter() const;
    [[nodiscard]] const Eigen::Tensor<Scalar,3> &get_M_bare() const;
    [[nodiscard]] const Eigen::Tensor<Scalar,3> &get_M()  const;
    [[nodiscard]] const Eigen::Tensor<Scalar,1> &get_L()  const;
    [[nodiscard]] const Eigen::Tensor<Scalar,1> &get_LC() const;
    [[nodiscard]] Eigen::Tensor<Scalar,3> &get_M_bare();
    [[nodiscard]] Eigen::Tensor<Scalar,3> &get_M();
    [[nodiscard]] Eigen::Tensor<Scalar,1> &get_L();
    [[nodiscard]] Eigen::Tensor<Scalar,1> &get_LC();
    [[nodiscard]] std::tuple<long,long,long>  get_dims() const;
    [[nodiscard]] long get_spin_dim() const;
    [[nodiscard]] long get_chiL()     const;
    [[nodiscard]] long get_chiR()     const;
    [[nodiscard]] size_t get_position() const;
    void set_position(size_t position_);
    void set_mps(const Eigen::Tensor<Scalar,3> &M_, const Eigen::Tensor<Scalar,1> &L_);

    void set_M(const Eigen::Tensor<Scalar,3> &M_);
    void set_L(const Eigen::Tensor<Scalar,1> &L_);
    void set_LC(const Eigen::Tensor<Scalar,1> &LC_);
    void unset_LC(){LC.reset();MC.reset();}
    void apply_mpo(const Eigen::Tensor<Scalar,4> & mpo);
    void apply_mpo(const Eigen::Tensor<Scalar,2> & mpo);

    //    [[nodiscard]] const Eigen::Tensor<Scalar,3> &get_G() const;
    //    [[nodiscard]] Eigen::Tensor<Scalar,3> &get_G();
    //    [[nodiscard]] Eigen::Tensor<Scalar,3> get_A()  const;
    //    [[nodiscard]] Eigen::Tensor<Scalar,3> get_B()  const;
    //    void set_mps_vidal(const Eigen::Tensor<Scalar,3> &G_, const Eigen::Tensor<Scalar,1> &L_);
    //    void set_mps_vidal(const Eigen::Tensor<Scalar,1> &L_,const Eigen::Tensor<Scalar,3> &G_);



};


#endif //DMRG_CLASS_MPS_SITE_H
