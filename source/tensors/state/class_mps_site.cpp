//
// Created by david on 2019-07-06.
//

#include "class_mps_site.h"
#include <general/nmspc_tensor_extra.h>
#include <io/nmspc_logger.h>

using Scalar = class_mps_site::Scalar;

class_mps_site::class_mps_site() = default;
class_mps_site::class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error)
    : M(M_), L(L_), position(pos), truncation_error(error){}



// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_mps_site::~class_mps_site() = default;                                              // default dtor
class_mps_site::class_mps_site(class_mps_site &&other)  noexcept = default;               // default move ctor
class_mps_site &class_mps_site::operator=(class_mps_site &&other) noexcept = default;     // default move assign
class_mps_site::class_mps_site(const class_mps_site &other) = default;
class_mps_site &class_mps_site::operator=(const class_mps_site &other) = default;








bool class_mps_site::isCenter() const {
    if(LC.has_value()) {
        if(LC.value().dimension(0) != M.dimension(2))
            throw std::runtime_error(fmt::format("M and LC dim mismatch: {} != {}", M.dimension(2), LC.value().dimension(0)));
    }
    return LC.has_value();
}





bool class_mps_site::is_real() const { return Textra::isReal(M, "M"); }
bool class_mps_site::has_nan() const { return Textra::hasNaN(M, "M"); }
void class_mps_site::assert_validity() const {
    if(Textra::hasNaN(M, "M")) throw std::runtime_error("MPS at position " + std::to_string(get_position()) + " has NAN's");
}

const Eigen::Tensor<Scalar, 3> &class_mps_site::get_M_bare() const { return std::as_const(M); }

const Eigen::Tensor<Scalar, 3> &class_mps_site::get_M() const {
    if(isCenter()) {
        if(MC) {
            return std::as_const(MC.value());
        } else {
            MC = M.contract(Textra::asDiagonal(get_LC()), Textra::idx({2}, {0}));
            return std::as_const(MC.value());
        }
    } else {
        return std::as_const(M);
    }
}

const Eigen::Tensor<Scalar, 1> &class_mps_site::get_L() const { return std::as_const(L); }
const Eigen::Tensor<Scalar, 1> &class_mps_site::get_LC() const {
    if(isCenter())
        return std::as_const(LC.value());
    else
        throw std::runtime_error("Site at position " + std::to_string(get_position()) + " is not a center");
}

Eigen::Tensor<Scalar, 3> &class_mps_site::get_M_bare() {
    return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M_bare());
}

Eigen::Tensor<Scalar, 3> &class_mps_site::get_M() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M()); }
Eigen::Tensor<Scalar, 1> &class_mps_site::get_L() { return L; }
Eigen::Tensor<Scalar, 1> &class_mps_site::get_LC() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_LC()); }

std::tuple<long, long, long> class_mps_site::get_dims() const { return {spin_dim(), get_chiL(), get_chiR()}; }
long                         class_mps_site::spin_dim() const { return M.dimension(0); }
long                         class_mps_site::get_chiL() const { return M.dimension(1); }
long                         class_mps_site::get_chiR() const { return M.dimension(2); }

void class_mps_site::set_position(const size_t position_) {
    position = position_;
    MC.reset();
}
size_t class_mps_site::get_position() const {
    if(position) {
        return position.value();
    } else {
        throw std::runtime_error("Position hasn't been set on mps site.");
    }
}

void class_mps_site::set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error) {
    set_M(M_);
    set_L(L_);
    set_truncation_error(error);
}

void class_mps_site::set_M(const Eigen::Tensor<Scalar, 3> &M_) {
    if(position) {
        M = M_;
        MC.reset();
    } else
        throw std::runtime_error("Can't set M: Position hasn't been set yet");
}
void class_mps_site::set_L(const Eigen::Tensor<Scalar, 1> &L_, double error) {
    if(position) {
        L                = L_;
        truncation_error = error;
    } else
        throw std::runtime_error("Can't set L: Position hasn't been set yet");
}
void class_mps_site::set_LC(const Eigen::Tensor<Scalar, 1> &LC_, double error) {
    if(position) {
        LC = LC_;
        MC.reset();
        truncation_error_LC = error;
    } else
        throw std::runtime_error("Can't set LC: Position hasn't been set yet");
}

void class_mps_site::set_truncation_error(double error) { truncation_error = error; }
void class_mps_site::set_truncation_error_LC(double error) { truncation_error_LC = error; }
double class_mps_site::get_truncation_error(){return truncation_error;}
double class_mps_site::get_truncation_error_LC(){return truncation_error_LC;}


void class_mps_site::apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo) {
    long mpoDimL = mpo.dimension(0);
    long mpoDimR = mpo.dimension(1);
    if(mpoDimL != mpoDimR) throw std::logic_error("Can't apply mpo's with different L/R dims: not implemented yet");

    if(isCenter()) {
        Eigen::Tensor<Scalar, 1> LC_temp = get_LC().broadcast(Textra::array1{mpoDimR});
        LC                               = LC_temp;
        MC.reset();
    }

    Eigen::Tensor<Scalar, 3> M_temp = M.contract(mpo, Textra::idx({0}, {2}))
                                          .shuffle(Textra::array5{4, 0, 2, 1, 3})
                                          .reshape(Textra::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    Eigen::Tensor<Scalar, 1> L_temp = L.broadcast(Textra::array1{mpoDimL});

    L = L_temp;
    M = M_temp;
}
void class_mps_site::apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo) {
    Eigen::Tensor<Scalar, 3> M_temp = mpo.contract(M, Textra::idx({0}, {0}));
    M                               = M_temp;
    MC.reset();
}

// void class_mps_site::set_mps_vidal  (const Eigen::Tensor<Scalar,3> &G_, const Eigen::Tensor<Scalar,1> &L_){set_G(G_); set_L(L_);}

// Eigen::Tensor<Scalar,3> class_mps_site::get_A()  const  {
//    return M;
//    return Textra::asDiagonal(L).contract(G, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
//}
// Eigen::Tensor<Scalar,3> class_mps_site::get_B()  const  {
//    return M;
//    return G.contract(Textra::asDiagonal(L), Textra::idx({2},{0}));
//}

// const Eigen::Tensor<Scalar,3> & class_mps_site::get_G() const {return std::as_const(G);}
//      Eigen::Tensor<Scalar,3> & class_mps_site::get_G()       {return G;}
