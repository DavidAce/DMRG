//
// Created by david on 2019-07-06.
//

#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
// textra must appear first
#include "class_mps_site.h"
#include <tools/common/fmt.h>

using Scalar = class_mps_site::Scalar;

class_mps_site::class_mps_site() = default;
class_mps_site::class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error)
    : M(M_), L(L_), position(pos), truncation_error(error) {}

class_mps_site::class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, std::optional<Eigen::Tensor<Scalar, 1>> L_, size_t pos, double error)
    : M(M_), L(std::move(L_)), position(pos), truncation_error(error) {}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_mps_site::~class_mps_site()                                      = default; // default dtor
class_mps_site::class_mps_site(class_mps_site &&other)                 = default; // default move ctor
class_mps_site &class_mps_site::operator=(class_mps_site &&other)      = default; // default move assign
class_mps_site::class_mps_site(const class_mps_site &other)            = default;
class_mps_site &class_mps_site::operator=(const class_mps_site &other) = default;

bool class_mps_site::isCenter() const {
    return LC.has_value();
}

bool class_mps_site::has_L() const { return L.has_value(); }
bool class_mps_site::has_M() const { return M.has_value(); }
bool class_mps_site::has_LC() const { return LC.has_value(); }

Eigen::DSizes<long, 3> class_mps_site::dimensions() const { return Eigen::DSizes<long, 3>{spin_dim(), get_chiL(), get_chiR()}; }

bool class_mps_site::is_real() const { return Textra::isReal(get_M_bare(), "M_bare") and Textra::isReal(get_L(), "L"); }
bool class_mps_site::has_nan() const { return Textra::hasNaN(get_M_bare(), "M_bare") and Textra::hasNaN(get_L(), "L"); }
void class_mps_site::assert_validity() const {
    if(has_nan()) throw std::runtime_error(fmt::format("class_mps_site::assert_validity(): MPS (M or L) at position {} has NaN's", get_position()));
}

const Eigen::Tensor<Scalar, 3> &class_mps_site::get_M_bare() const {
    if(M) {
        if(M.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_M_bare(): M has size 0 at position {}", get_position()));
        return M.value();
    } else
        throw std::runtime_error(fmt::format("class_mps_site::get_M_bare(): M has not been set at position {}", get_position()));
}
const Eigen::Tensor<Scalar, 3> &class_mps_site::get_M() const {
    if(isCenter()) {
        if(LC.value().dimension(0) != get_M_bare().dimension(2))
            throw std::runtime_error(
                fmt::format("class_mps_site::get_M(): M and LC dim mismatch: {} != {}", get_M_bare().dimension(2), LC.value().dimension(0)));
        if(MC){
            if(MC.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_M(): MC has size 0 at position {}", get_position()));
            return MC.value();
        } else {
            MC = get_M_bare().contract(Textra::asDiagonal(get_LC()), Textra::idx({2}, {0}));
            if(MC.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_M(): built MC with size 0 at position {}", get_position()));
            return MC.value();
        }
    } else
        return get_M_bare();
}

const Eigen::Tensor<Scalar, 1> &class_mps_site::get_L() const {
    if(L) {
        if(L.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_L(): L has size 0 at position {}", get_position()));
        return L.value();
    } else
        throw std::runtime_error(fmt::format("class_mps_site::get_L(): L has not been set at position {}", get_position()));
}
const Eigen::Tensor<Scalar, 1> &class_mps_site::get_LC() const {
    if(isCenter()) {
        if(LC.value().dimension(0) != get_M_bare().dimension(2))
            throw std::runtime_error(
                fmt::format("class_mps_site::get_LC(): M and LC dim mismatch: {} != {}", get_M_bare().dimension(2), LC.value().dimension(0)));
        if(LC.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_LC(): LC has size 0 at position {}", get_position()));
        return LC.value();
    } else
        throw std::runtime_error(fmt::format("class_mps_site::get_LC(): Site at position {} is not a center", get_position()));
}

Eigen::Tensor<Scalar, 3> &   class_mps_site::get_M_bare() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M_bare()); }
Eigen::Tensor<Scalar, 3> &   class_mps_site::get_M() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M()); }
Eigen::Tensor<Scalar, 1> &   class_mps_site::get_L() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_L()); }
Eigen::Tensor<Scalar, 1> &   class_mps_site::get_LC() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_LC()); }
std::tuple<long, long, long> class_mps_site::get_dims() const { return {spin_dim(), get_chiL(), get_chiR()}; }
long                         class_mps_site::spin_dim() const { return get_M_bare().dimension(0); }
long                         class_mps_site::get_chiL() const { return get_M_bare().dimension(1); }
long                         class_mps_site::get_chiR() const { return get_M_bare().dimension(2); }

void class_mps_site::set_position(const size_t position_) {
    position = position_;
    MC.reset();
}
size_t class_mps_site::get_position() const {
    if(position) {
        return position.value();
    } else {
        throw std::runtime_error("class_mps_site::get_position(): Position hasn't been set on mps site.");
    }
}

void class_mps_site::set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error) {
    // M has to be a "bare" matrix, i.e. not an MC which would include LC.
    set_M(M_);
    set_L(L_);
    set_truncation_error(error);
}

void class_mps_site::set_M(const Eigen::Tensor<Scalar, 3> &M_) {
    // M has to be a "bare" matrix, i.e. not an MC which would include LC.
    if(position) {
        M = M_;
        MC.reset();
    } else
        throw std::runtime_error("class_mps_site::set_M(const Eigen::Tensor<Scalar, 3> &): Can't set M: Position hasn't been set yet");
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

void   class_mps_site::set_truncation_error(double error) { truncation_error = error; }
void   class_mps_site::set_truncation_error_LC(double error) { truncation_error_LC = error; }
double class_mps_site::get_truncation_error() const { return truncation_error; }
double class_mps_site::get_truncation_error_LC() const { return truncation_error_LC; }

void class_mps_site::merge_mps(const class_mps_site &other) {
    // This operation is done when merging mps after an svd split, for instance
    if(get_position() != other.get_position())
        throw std::runtime_error(fmt::format("class_mps_site::merge_mps(const class_mps_site &): position mismatch: this pos {} != other pos {}",
                                             get_position(), other.get_position()));

    if(not other.has_M()) throw std::runtime_error("class_mps_site::merge_mps(const class_mps_site &): Got mps site with undefined M");

    // We have to copy the bare "M", i.e. not MC, which would include LC.
    set_M(other.get_M_bare());

    // We should only copy the L matrix if it exists
    // If it does not exist it may just not have been set.
    // For instance, when splitting a multisite tensor into its contituent mps sites,
    // the edge mps sites have empty L matrices on them because these are not updated.
    // This is intended, but should be checked for.
    // Note that if there is no L on the other mps, there has to exist a previous one
    // on this mps -- we can't leave here without having defined L.
    if(other.has_L())
        set_L(other.get_L(), other.get_truncation_error());
    else if(not has_L())
        throw std::runtime_error(fmt::format(
            "class_mps_site::merge_mps(const class_mps_site &): Got mps site with undefined L, and no L is defined on this mps either, at position {}",
            get_position()));

    // If the other is a center it should be fine to copy its contents
    // without checking its size. However for debugging purposes, since
    // that would be unexpected, just give a warning.
    if(other.isCenter()) {
        set_LC(other.get_LC(), other.get_truncation_error_LC());
        if(get_M_bare().dimension(2) != get_LC().dimension(0))
            throw std::runtime_error(
                fmt::format("Got center of wrong dimension. M dims {} are not compatible with LC dims {}", get_M_bare().dimensions(), get_LC().dimensions()));
    }
}

void class_mps_site::apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo) {
    long mpoDimL = mpo.dimension(0);
    long mpoDimR = mpo.dimension(1);
    if(mpoDimL != mpoDimR) throw std::logic_error("Can't apply mpo's with different L/R dims: not implemented yet");

    if(isCenter()) {
        Eigen::Tensor<Scalar, 1> LC_temp = get_LC().broadcast(Textra::array1{mpoDimR});
        set_LC(LC_temp);
    }
    Eigen::Tensor<Scalar, 3> M_bare_temp(Textra::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    M_bare_temp.device(Textra::omp::getDevice()) = get_M_bare()
                                               .contract(mpo, Textra::idx({0}, {2}))
                                               .shuffle(Textra::array5{4, 0, 2, 1, 3})
                                               .reshape(Textra::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    Eigen::Tensor<Scalar, 1> L_temp = get_L().broadcast(Textra::array1{mpoDimL});

    set_L(L_temp);
    set_M(M_bare_temp);
}

void class_mps_site::apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo) {
    Eigen::Tensor<Scalar, 3> M_bare_temp(mpo.dimension(1), get_chiL(), get_chiR());
    M_bare_temp.device(Textra::omp::getDevice()) = mpo.contract(get_M_bare(), Textra::idx({0}, {0}));
    set_M(M_bare_temp);
}
