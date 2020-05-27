//
// Created by david on 2017-11-13.
//

#include <config/nmspc_settings.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_mps_site.h>
#include <tools/common/log.h>
using namespace std;
using namespace Textra;

using Scalar = class_mps_2site::Scalar;

class_mps_2site::class_mps_2site()
    : MPS_A(std::make_unique<class_mps_site>()),
      MPS_B(std::make_unique<class_mps_site>())
{
    tools::log->trace("Constructing 2site mps");
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
class_mps_2site::~class_mps_2site() = default;                                                    // default dtor
class_mps_2site::class_mps_2site(class_mps_2site &&other)  noexcept = default;               // default move ctor
class_mps_2site &class_mps_2site::operator=(class_mps_2site &&other) noexcept = default;     // default move assign

class_mps_2site::class_mps_2site(const class_mps_2site &other):
    MPS_A(std::make_unique<class_mps_site>(*other.MPS_A)),
    MPS_B(std::make_unique<class_mps_site>(*other.MPS_B)),
    swapped(other.swapped)
{}

class_mps_2site &class_mps_2site::operator=(const class_mps_2site &other) {
    // check for self-assignment
    if(this != &other) {
        MPS_A = std::make_unique<class_mps_site>(*other.MPS_A);
        MPS_B = std::make_unique<class_mps_site>(*other.MPS_B);
        swapped = other.swapped;
    }
    return *this;
}


void class_mps_2site::initialize(long spin_dim) {
    tools::log->trace("Initializing 2site mps");
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

void class_mps_2site::assert_validity() const {
    MPS_A->assert_validity();
    MPS_B->assert_validity();
}
bool class_mps_2site::is_real() const { return MPS_A->is_real() and MPS_B->is_real(); }
bool class_mps_2site::has_nan() const { return MPS_A->has_nan() or MPS_B->has_nan(); }
long class_mps_2site::chiC() const { return MPS_A->get_LC().dimension(0); }
long class_mps_2site::chiA() const { return MPS_A->get_L().dimension(0); }
long class_mps_2site::chiB() const { return MPS_B->get_L().dimension(0); }
long class_mps_2site::spin_dim_A() const { return MPS_A->spin_dim(); }
long class_mps_2site::spin_dim_B() const { return MPS_B->spin_dim(); }

const Eigen::Tensor<Scalar, 3> &class_mps_2site::A_bare() const { return MPS_A->get_M_bare(); }
const Eigen::Tensor<Scalar, 3> &class_mps_2site::A() const { return MPS_A->get_M(); }
const Eigen::Tensor<Scalar, 3> &class_mps_2site::B() const { return MPS_B->get_M(); }
Eigen::Tensor<Scalar, 2>        class_mps_2site::LC() const { return Textra::asDiagonal(MPS_A->get_LC()); }

Eigen::Tensor<Scalar, 3> class_mps_2site::GA() const {
    return MPS_A->get_M_bare().contract(Textra::asDiagonalInversed(MPS_A->get_L()), Textra::idx({1}, {1})).shuffle(Textra::array3{0, 2, 1});
}
Eigen::Tensor<Scalar, 3> class_mps_2site::GB() const { return MPS_B->get_M_bare().contract(Textra::asDiagonalInversed(MPS_B->get_L()), Textra::idx({2}, {0})); }
Eigen::Tensor<Scalar, 2> class_mps_2site::LA() const { return Textra::asDiagonal(MPS_A->get_L()); }
Eigen::Tensor<Scalar, 2> class_mps_2site::LB() const { return Textra::asDiagonal(MPS_B->get_L()); }

void class_mps_2site::set_mps(const class_mps_site & mpsA, const class_mps_site & mpsB){
    if(not mpsA.isCenter())
        throw std::runtime_error("Given mps for site A is not a center");
    MPS_A = std::make_unique<class_mps_site>(mpsA);
    MPS_B = std::make_unique<class_mps_site>(mpsB);
}



void class_mps_2site::set_mps(const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB) {
    MPS_A->set_M(MA);
    MPS_A->set_LC(LC);
    MPS_B->set_M(MB);
}

void class_mps_2site::set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC,
                              const Eigen::Tensor<Scalar, 3> &MB, const Eigen::Tensor<Scalar, 1> &LB) {
    MPS_A->set_mps(MA, LA);
    MPS_A->set_LC(LC);
    MPS_B->set_mps(MB, LB);
}

Eigen::DSizes<long, 3> class_mps_2site::dimensions() const { return Eigen::DSizes<long, 3>{spin_dim_A() * spin_dim_B(), chiA(), chiB()}; }

void class_mps_2site::swap_AB() {
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

const class_mps_site &class_mps_2site::get_mps_siteA() const { return *MPS_A; }
const class_mps_site &class_mps_2site::get_mps_siteB() const { return *MPS_B; }
class_mps_site &      class_mps_2site::get_mps_siteA() { return *MPS_A; }
class_mps_site &      class_mps_2site::get_mps_siteB() { return *MPS_B; }

void class_mps_2site::set_positions(size_t posA, size_t posB) {
    MPS_A->set_position(posA);
    MPS_B->set_position(posB);
}
void class_mps_2site::set_positionA(size_t pos) { MPS_A->set_position(pos); }
void class_mps_2site::set_positionB(size_t pos) { MPS_B->set_position(pos); }

std::pair<size_t, size_t> class_mps_2site::get_positions() const { return std::make_pair(MPS_A->get_position(), MPS_B->get_position()); }
size_t                    class_mps_2site::get_positionA() const { return MPS_A->get_position(); }
size_t                    class_mps_2site::get_positionB() const { return MPS_B->get_position(); }

Eigen::Tensor<class_mps_2site::Scalar, 3> class_mps_2site::get_2site_tensor(Scalar norm) const
/*!
 * Returns a two-site tensor
     @verbatim
        1--[ LA_diag ]--[ GA ]--[ LC_diag ]-- [ GB ] -- [ LB_diag ]--3
                     |                 |
                     0                 2
       which in our notation is simply A * B  (because A = LA_diag * GA * LC_diag and B = GB * LB_diag),
       becomes

        1--[ 2site_tensor ]--2
                  |
                  0

     @endverbatim
 */
{
    long dim0 = MPS_A->spin_dim() * MPS_B->spin_dim();
    long dim1 = MPS_A->get_chiL();
    long dim2 = MPS_B->get_chiR();
    return A().contract(B(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dim0, dim1, dim2}) / norm;
}
