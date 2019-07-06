//
// Created by david on 2019-07-06.
//

#include "class_vidal_site.h"

using Scalar = class_vidal_site::Scalar;
bool class_vidal_site::isReal()const{return Textra::isReal(G,"G");}
const Eigen::Tensor<Scalar,3> & class_vidal_site::get_G() const {return std::as_const(G);}
const Eigen::Tensor<Scalar,1> & class_vidal_site::get_L() const {return std::as_const(L);}
      Eigen::Tensor<Scalar,3> & class_vidal_site::get_G()       {return G;}
      Eigen::Tensor<Scalar,1> & class_vidal_site::get_L()       {return L;}



Eigen::Tensor<Scalar,3> class_vidal_site::get_A()  const  {return Textra::asDiagonal(L).contract(G, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});}
Eigen::Tensor<Scalar,3> class_vidal_site::get_B()  const  {return G.contract(Textra::asDiagonal(L), Textra::idx({2},{0}));}

std::tuple<long,long,long> class_vidal_site::get_dims() const {return {get_spin_dim(),get_chiL(),get_chiR()};}
long class_vidal_site::get_spin_dim() const {return G.dimension(0);}
long class_vidal_site::get_chiL()     const {return G.dimension(1);}
long class_vidal_site::get_chiR()     const {return G.dimension(2);}



void   class_vidal_site::set_position(const size_t position_){position = position_;}
size_t class_vidal_site::get_position() const {
    if(position) {return position.value();}
    else{throw std::runtime_error("Position hasn't been set on vidal site.");}
}



void class_vidal_site::set_mps  (const Eigen::Tensor<Scalar,3> &G_, const Eigen::Tensor<Scalar,1> &L_){set_G(G_); set_L(L_);}
void class_vidal_site::set_L    (const Eigen::Tensor<Scalar,1> &L_){if(position) L=L_; else throw std::runtime_error("Can't set L: Position hasn't been set yet");}
void class_vidal_site::set_G    (const Eigen::Tensor<Scalar,3> &G_){if(position) G=G_; else throw std::runtime_error("Can't set G: Position hasn't been set yet");}

