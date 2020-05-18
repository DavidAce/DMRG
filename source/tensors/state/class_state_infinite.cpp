//
// Created by david on 7/22/17.
//

#include <config/nmspc_settings.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/common/svd.h>
#include <tools/common/views.h>
#include <tools/infinite/measure.h>

using Scalar = class_state_infinite::Scalar;

class_state_infinite::class_state_infinite() : MPS_A(std::make_unique<class_mps_site>()), MPS_B(std::make_unique<class_mps_site>()) {
    tools::log->trace("Constructing state");
}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_state_infinite::~class_state_infinite()                                     = default;            // default dtor
class_state_infinite::class_state_infinite(class_state_infinite &&other) noexcept = default;            // default move ctor
class_state_infinite &class_state_infinite::operator=(class_state_infinite &&other) noexcept = default; // default move assign

/* clang-format off */
class_state_infinite::class_state_infinite(const class_state_infinite &other):
MPS_A(std::make_unique<class_mps_site>(*other.MPS_A)),
MPS_B(std::make_unique<class_mps_site>(*other.MPS_B)),
swapped(other.swapped),
chi_lim(other.chi_lim),
chi_max(other.chi_max),
cache(other.cache),
measurements(other.measurements),
lowest_recorded_variance(other.lowest_recorded_variance)
{}
/* clang-format on */

class_state_infinite &class_state_infinite::operator=(const class_state_infinite &other) {
    // check for self-assignment
    if(this != &other) {
        MPS_A                    = std::make_unique<class_mps_site>(*other.MPS_A);
        MPS_B                    = std::make_unique<class_mps_site>(*other.MPS_B);
        swapped                  = other.swapped;
        chi_lim                  = other.chi_lim;
        chi_max                  = other.chi_max;
        cache                    = other.cache;
        measurements             = other.measurements;
        lowest_recorded_variance = other.lowest_recorded_variance;
    }
    return *this;
}

void class_state_infinite::initialize(ModelType model_type) {
    tools::log->trace("Initializing state");
    long spin_dim = 2;
    switch(model_type) {
        case ModelType::ising_tf_rf: spin_dim = settings::model::ising_tf_rf::spin_dim; break;
        case ModelType::ising_sdual: spin_dim = settings::model::ising_sdual::spin_dim; break;
        default: spin_dim = 2;
    }
    Eigen::Tensor<Scalar, 3> M(spin_dim, 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);
    // Default is a product state, spins pointing up in z.
    M.setZero();
    M(0, 0, 0) = 1;
    L.setConstant(1.0);
    MPS_A = std::make_unique<class_mps_site>(M, L, 0);
    MPS_B = std::make_unique<class_mps_site>(M, L, 1);
    MPS_A->set_LC(L);
}

std::pair<size_t, size_t> class_state_infinite::get_positions() { return std::make_pair(MPS_A->get_position(), MPS_B->get_position()); }
size_t                    class_state_infinite::get_positionA() { return MPS_A->get_position(); }
size_t                    class_state_infinite::get_positionB() { return MPS_B->get_position(); }

long class_state_infinite::get_chi() const { return MPS_A->get_LC().dimension(0); }
long class_state_infinite::get_chiA() const { return MPS_A->get_L().dimension(0); }
long class_state_infinite::get_chiB() const { return MPS_B->get_L().dimension(0); }

long class_state_infinite::get_chi_lim() const {
    // Should get the the current limit on allowed bond dimension
    return chi_lim.value();
}
long class_state_infinite::get_chi_max() const {
    // Should get the the current limit on allowed bond dimension for the duration of the simulation
    return chi_max.value();
}

long class_state_infinite::get_spin_dimA() const { return MPS_A->spin_dim(); }
long class_state_infinite::get_spin_dimB() const { return MPS_B->spin_dim(); }

void class_state_infinite::set_chi_lim(long chi_lim_) {
    // Should set the the current limit on allowed bond dimension
    if(chi_lim_ == 0) throw std::runtime_error("Can't set chi limit to zero!");
    chi_lim = chi_lim_;
}

void class_state_infinite::set_chi_max(long chi_max_) {
    // Should set the the current limit on allowed bond dimension for the duration of the simulation
    if(chi_max_ == 0) throw std::runtime_error("Can't set chi max to zero!");
    chi_max = chi_max_;
}

double class_state_infinite::get_truncation_error() const {
    // Should get the the current limit on allowed bond dimension for the duration of the simulation
    return tools::infinite::measure::truncation_error(*this);
}

Eigen::DSizes<long, 3> class_state_infinite::dimensions() const { return Eigen::DSizes<long, 3>{get_spin_dimA() * get_spin_dimB(), get_chiA(), get_chiB()}; }

const class_mps_site &class_state_infinite::get_mps_siteA() const { return *MPS_A; }
const class_mps_site &class_state_infinite::get_mps_siteB() const { return *MPS_B; }
class_mps_site &      class_state_infinite::get_mps_siteA() { return *MPS_A; }
class_mps_site &      class_state_infinite::get_mps_siteB() { return *MPS_B; }

const Eigen::Tensor<Scalar, 3> &class_state_infinite::A_bare() const { return MPS_A->get_M_bare(); }
const Eigen::Tensor<Scalar, 3> &class_state_infinite::A() const { return MPS_A->get_M(); }
const Eigen::Tensor<Scalar, 3> &class_state_infinite::B() const { return MPS_B->get_M(); }
Eigen::Tensor<Scalar, 2>        class_state_infinite::LC() const { return Textra::asDiagonal(MPS_A->get_LC()); }
Eigen::Tensor<Scalar, 3>        class_state_infinite::GA() const {
    return MPS_A->get_M_bare().contract(Textra::asDiagonalInversed(MPS_A->get_L()), Textra::idx({1}, {1})).shuffle(Textra::array3{0, 2, 1});
}
Eigen::Tensor<Scalar, 3> class_state_infinite::GB() const {
    return MPS_B->get_M_bare().contract(Textra::asDiagonalInversed(MPS_B->get_L()), Textra::idx({2}, {0}));
}
Eigen::Tensor<Scalar, 2> class_state_infinite::LA() const { return Textra::asDiagonal(MPS_A->get_L()); }
Eigen::Tensor<Scalar, 2> class_state_infinite::LB() const { return Textra::asDiagonal(MPS_B->get_L()); }

const Eigen::Tensor<Scalar, 3> &class_state_infinite::get_2site_tensor(Scalar norm) const {
    /*!
 * Returns a two-site tensor
     @verbatim
        1--[ LA ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
       which in our notation is simply A * B  (because A = LA * GA * LC and B = GB * LB),
       becomes

        1--[ 2site_tensor ]--2
                  |
                  0

     @endverbatim
 */

    if(cache.twosite_tensor) return cache.twosite_tensor.value();
    long dim0            = MPS_A->spin_dim() * MPS_B->spin_dim();
    long dim1            = MPS_A->get_chiL();
    long dim2            = MPS_B->get_chiR();
    cache.twosite_tensor = A().contract(B(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dim0, dim1, dim2}) / norm;
    return cache.twosite_tensor.value();
}

void class_state_infinite::assert_validity() const {
    MPS_A->assert_validity();
    MPS_B->assert_validity();
}
bool class_state_infinite::is_real() const { return MPS_A->is_real() and MPS_B->is_real(); }
bool class_state_infinite::has_nan() const { return MPS_A->has_nan() or MPS_B->has_nan(); }
//
// template<typename T>
// Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_matrix() const {
//    Eigen::Tensor<T, 5> tempL;
//    Eigen::Tensor<T, 5> tempR;
//    if constexpr(std::is_same<T, double>::value) {
//        if(not Lblock->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Lblock when building H_local");
//        }
//        if(not Rblock->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Rblock when building H_local");
//        }
//        if(not HA->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO A when building H_local");
//        }
//        if(not HB->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO B when building H_local");
//        }
//        tempL = Lblock->block.contract(HA->MPO(), Textra::idx({2}, {0})).real().shuffle(Textra::array5{4, 1, 3, 0, 2}).real();
//        tempR = Rblock->block.contract(HB->MPO(), Textra::idx({2}, {1})).real().shuffle(Textra::array5{4, 1, 3, 0, 2}).real();
//    } else {
//        tempL = Lblock->block.contract(HA->MPO(), Textra::idx({2}, {0})).shuffle(Textra::array5{4, 1, 3, 0, 2});
//        tempR = Rblock->block.contract(HB->MPO(), Textra::idx({2}, {1})).shuffle(Textra::array5{4, 1, 3, 0, 2});
//    }
//    long                shape   = mps_sites->chiA() * mps_sites->spindim() * mps_sites->chiB() * mps_sites->spin_dim_A();
//    Eigen::Tensor<T, 8> H_local = tempL.contract(tempR, Textra::idx({4}, {4})).shuffle(Textra::array8{0, 1, 4, 5, 2, 3, 6, 7});
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(H_local.data(), shape, shape);
//}

// template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>               class_state_infinite::get_H_local_matrix<double>() const;
// template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_matrix<std::complex<double>>() const;
//
// template<typename T>
// Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_sq_matrix() const {
//    Eigen::Tensor<T, 6> tempL;
//    Eigen::Tensor<T, 6> tempR;
//    if constexpr(std::is_same<T, double>::value) {
//        if(not Lblock2->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Lblock2 when building H_local_sq");
//        }
//        if(not Rblock2->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Rblock2 when building H_local_sq");
//        }
//        if(not HA->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO A when building H_local_sq");
//        }
//        if(not HB->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO B when building H_local_sq");
//        }
//        tempL = Lblock2->block.contract(HA->MPO(), Textra::idx({2}, {0}))
//                    .contract(HA->MPO(), Textra::idx({2, 5}, {0, 2}))
//                    .real()
//                    .shuffle(Textra::array6{5, 1, 3, 0, 2, 4});
//        tempR = Rblock2->block.contract(HB->MPO(), Textra::idx({2}, {1}))
//                    .contract(HB->MPO(), Textra::idx({2, 5}, {1, 2}))
//                    .real()
//                    .shuffle(Textra::array6{5, 1, 3, 0, 2, 4});
//    } else {
//        tempL = Lblock2->block.contract(HA->MPO(), Textra::idx({2}, {0}))
//                    .contract(HA->MPO(), Textra::idx({2, 5}, {0, 2}))
//                    .shuffle(Textra::array6{5, 1, 3, 0, 2, 4});
//        tempR = Rblock2->block.contract(HB->MPO(), Textra::idx({2}, {1}))
//                    .contract(HB->MPO(), Textra::idx({2, 5}, {1, 2}))
//                    .shuffle(Textra::array6{5, 1, 3, 0, 2, 4});
//    }
//
//    long                shape   = mps_sites->chiA() * mps_sites->spindim() * mps_sites->chiB() * mps_sites->spin_dim_A();
//    Eigen::Tensor<T, 8> H_local = tempL.contract(tempR, Textra::idx({4, 5}, {4, 5})).shuffle(Textra::array8{0, 1, 4, 5, 2, 3, 6, 7});
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(H_local.data(), shape, shape);
//}

// template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>               class_state_infinite::get_H_local_sq_matrix<double>() const;
// template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_sq_matrix<std::complex<double>>() const;

// void class_state_infinite::enlarge_environment(int direction) {
//    assert_positions();
//    auto position_A_new = mps_sites->MPS_A->get_position() + 1;
//    *Lblock             = Lblock->enlarge(*mps_sites->MPS_A, *HA);
//    *Rblock             = Rblock->enlarge(*mps_sites->MPS_B, *HB);
//    *Lblock2            = Lblock2->enlarge(*mps_sites->MPS_A, *HA);
//    *Rblock2            = Rblock2->enlarge(*mps_sites->MPS_B, *HB);
//    set_positions(position_A_new);
//    if(direction != 0) throw std::runtime_error("Ooops, direction != 0");
//    if (direction == 1){
//        *Lblock  = Lblock->enlarge(*mps_sites->MPS_A,  *HA);
//        *Lblock2 = Lblock2->enlarge(*mps_sites->MPS_A, *HA);
//        Lblock->set_position (HB->get_position());
//        Lblock2->set_position(HB->get_position());
//    }else if (direction == -1){
//        *Rblock  = Rblock->enlarge(*mps_sites->MPS_B,  *HB);
//        *Rblock2 = Rblock2->enlarge(*mps_sites->MPS_B, *HB);
//        Rblock->set_position (HA->get_position());
//        Rblock2->set_position(HA->get_position());
//    }else if(direction == 0){
//
//
//    }
//    clear_measurements();
//    assert_positions();
//}

void class_state_infinite::set_positions(size_t position) {
    MPS_A->set_position(position);
    MPS_B->set_position(position + 1);
}

void class_state_infinite::set_mps(const Eigen::Tensor<Scalar, 3> &twosite_tensor) {
    long   chi      = get_chi_lim();
    long   dA       = get_spin_dimA();
    long   dB       = get_spin_dimB();
    size_t posA     = get_positionA();
    size_t posB     = get_positionB();
    auto   mps_list = tools::common::svd::split_mps(twosite_tensor, {dA, dB}, {posA, posB}, posA, chi);
    set_mps(mps_list);
}

void class_state_infinite::set_mps(const std::list<class_mps_site> &mps_list) {
    if(mps_list.size() != 2) throw std::runtime_error("Expected 2 sites, got: " + std::to_string(mps_list.size()));
    const auto &mpsA = *std::next(mps_list.begin(), 0);
    const auto &mpsB = *std::next(mps_list.begin(), 1);
    set_mps(mpsA, mpsB);
    clear_cache();
}

void class_state_infinite::set_mps(const class_mps_site &mpsA, const class_mps_site &mpsB) {
    if(not mpsA.isCenter()) throw std::runtime_error("Given mps for site A is not a center");
    MPS_A = std::make_unique<class_mps_site>(mpsA);
    MPS_B = std::make_unique<class_mps_site>(mpsB);
    clear_cache();
}

void class_state_infinite::set_mps(const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB) {
    MPS_A->set_M(MA);
    MPS_A->set_LC(LC);
    MPS_B->set_M(MB);
    clear_cache();
}
void class_state_infinite::set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC,
                                   const Eigen::Tensor<Scalar, 3> &MB, const Eigen::Tensor<Scalar, 1> &LB) {
    MPS_A->set_mps(MA, LA);
    MPS_A->set_LC(LC);
    MPS_B->set_mps(MB, LB);
    clear_cache();
}

void class_state_infinite::clear_cache() const { cache = Cache(); }

void class_state_infinite::clear_measurements() const {
    measurements                              = state_measure_infinite();
    tools::common::views::components_computed = false;
}

void class_state_infinite::do_all_measurements() const { tools::infinite::measure::do_all_measurements(*this); }

void class_state_infinite::swap_AB() {
    tools::log->trace("Swapping AB");
    swapped = !swapped;
    // Store the positions
    auto position_left  = MPS_A->get_position();
    auto position_right = MPS_B->get_position();

    // Swap Gamma
    Eigen::Tensor<Scalar, 1> LC = MPS_A->get_LC();
    MPS_A->unset_LC();
    MPS_B->unset_LC();
    MPS_A.swap(MPS_B);
    MPS_A->set_LC(MPS_A->get_L());
    MPS_A->set_L(LC);
    MPS_B->set_L(LC);
    MPS_A->set_position(position_left);
    MPS_B->set_position(position_right);
}
