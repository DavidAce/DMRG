//
// Created by david on 2017-11-13.
//


#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>
#include <math/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <simulation/nmspc_settings.h>
#include <iomanip>

using namespace std;
using namespace Textra;

using Scalar = class_mps_2site::Scalar;

// Function definitions

class_mps_2site::class_mps_2site(){
    MPS_A = std::make_unique<class_mps_site>();
    MPS_B = std::make_unique<class_mps_site>();
}

class_mps_2site::class_mps_2site(const class_mps_2site &other)
    : spin_dimension(other.spin_dimension), swapped(other.swapped), MPS_A(copy_unique(other.MPS_A)), MPS_B(copy_unique(other.MPS_B))
{
}


bool class_mps_2site::isReal()const {return MPS_A->isReal() and MPS_B->isReal();}
long class_mps_2site::chiA () const {return MPS_A->get_L().dimension(0);}
long class_mps_2site::chiB () const {return MPS_B->get_L().dimension(0);}
long class_mps_2site::chiC () const {return MPS_A->get_LC().dimension(0);}
long class_mps_2site::spindim() const {return spin_dimension;}
const Eigen::Tensor<Scalar,3> & class_mps_2site::A() const {return MPS_A->get_M();}
const Eigen::Tensor<Scalar,3> & class_mps_2site::B() const {return MPS_B->get_M();}
Eigen::Tensor<Scalar,2> class_mps_2site::C() const {return Textra::asDiagonal(MPS_A->get_LC());}

Eigen::Tensor<Scalar,3> class_mps_2site::GA() const {
    return MPS_A->get_M_bare()
    .contract(Textra::asDiagonalInversed(MPS_A->get_L()), Textra::idx({1},{1}))
    .shuffle(Textra::array3{0,2,1});
}
Eigen::Tensor<Scalar,3> class_mps_2site::GB() const {
    return MPS_B->get_M()
    .contract(Textra::asDiagonalInversed(MPS_B->get_L()), Textra::idx({2},{0}));
}
Eigen::Tensor<Scalar,2> class_mps_2site::LA() const {return Textra::asDiagonal(MPS_A->get_L());}
Eigen::Tensor<Scalar,2> class_mps_2site::LB() const {return Textra::asDiagonal(MPS_B->get_L());}

void class_mps_2site::set_mps(const Eigen::Tensor<Scalar,1> &LA,
                              const Eigen::Tensor<Scalar,3> &U,
                              const Eigen::Tensor<Scalar,1> &LC,
                              const Eigen::Tensor<Scalar,3> &V_dagger,
                              const Eigen::Tensor<Scalar,1> &LB)
{
    MPS_A->set_mps(U, LA);
    MPS_A->set_LC(LC);
    MPS_B->set_mps(V_dagger, LB);
}

Eigen::DSizes<long,4> class_mps_2site::dimensions() const {return Eigen::DSizes<long,4>{spin_dimension,chiA(), spin_dimension,chiB()};}



void class_mps_2site::initialize(int spin_dim){
    spin_dimension = spin_dim;
    Eigen::Tensor<Scalar,3> UA(spin_dimension,1,1);
    Eigen::Tensor<Scalar,3> VB(spin_dimension,1,1);
    Eigen::Tensor<Scalar,1> LA(1);
    Eigen::Tensor<Scalar,1> LC(1);
    Eigen::Tensor<Scalar,1> LB(1);

    // Default is a product state, spins pointing up in z.
    UA.setZero();
    VB.setZero();
    LA.setConstant(1.0);
    LB.setConstant(1.0);
    LC.setConstant(1.0);

    UA(0,0,0) = 1;
    VB(0,0,0) = 1;

    set_mps(LA,UA,LC,VB,LB);
    swapped = false;
}


void class_mps_2site::swap_AB() {
    swapped = !swapped;

    //Swap Gamma
    MPS_A.swap(MPS_B);
    MPS_A->set_LC(MPS_A->get_L());
    MPS_A->set_L (MPS_B->get_LC());
    MPS_B->unset_LC(); // Delete LC from B.
//
//    tmp3 = MPS_A->get_G();
//    MPS_A->set_G(MPS_B->get_G());
//    MPS_B->set_G(tmp3);
//
//
//    tmp1 = LC;
//    LC = MPS_B->get_L();
//    MPS_A->set_L(tmp1);
//    MPS_B->set_L(tmp1);
}


Eigen::Tensor<class_mps_2site::Scalar,4> class_mps_2site::get_theta(Scalar norm) const
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LA ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return
            A().contract(B(), Textra::idx({2},{1})) / norm;
//    return
//            A().contract(C(), idx({2},{0}))
//               .contract(B(), idx({2},{1})) / norm;
}
