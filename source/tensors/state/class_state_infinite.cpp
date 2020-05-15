//
// Created by david on 7/22/17.
//

//#include <tensors/state/class_optimize_mps.h>
#include <config/nmspc_settings.h>
#include <iomanip>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <math/class_eigsolver.h>
#include <math/nmspc_math.h>
#include <tensors/model/class_mpo_factory.h>
#include <tensors/state/class_environment.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/common/views.h>
#include <tools/infinite/measure.h>

#define profile_optimization 0

using namespace std;
using namespace Textra;
using Scalar = class_state_infinite::Scalar;

class_state_infinite::class_state_infinite()
    : MPS(std::make_unique<class_mps_2site>())
//        Lblock2(std::make_unique<class_environment_var>("L",0)),
//        Rblock2(std::make_unique<class_environment_var>("R",1))
{
    tools::log->trace("Constructing class_state_infinite");
}

// We need to make a destructor manually for the enclosing class "class_state_infinite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
// class_state_infinite::~class_state_infinite()=default;

// we could the default move constructor
// class_state_infinite::class_state_infinite(class_state_infinite&&)  = default;
// we could use the default move assignment operator
// class_state_infinite & class_state_infinite::operator=(class_state_infinite&&) = default;

// ...but we would have to provide a deep copying assignment operator
// class_state_infinite & class_state_infinite::operator=(class_state_infinite const& source) {
//    HA = source.HA->clone();
//    HB = source.HB->clone();
//    return *this;
//}


long class_state_infinite::get_chi() const { return MPS->chiC(); }

long class_state_infinite::get_chi_lim() const {
    // Should get the the current limit on allowed bond dimension
    return chi_lim.value();
}
void class_state_infinite::set_chi_lim(long chi_lim_) {
    // Should set the the current limit on allowed bond dimension
    if(chi_lim_ == 0) throw std::runtime_error("Can't set chi limit to zero!");
    chi_lim = chi_lim_;
}

long class_state_infinite::get_chi_max() const {
    // Should get the the current limit on allowed bond dimension for the duration of the simulation
    return chi_max.value();
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

Eigen::DSizes<long, 3> class_state_infinite::dimensions() const { return MPS->dimensions(); }

Eigen::Tensor<Scalar, 3> class_state_infinite::get_mps() const {
    if(cache.mps) return cache.mps.value();
    cache.mps = MPS->get_mps();
    return cache.mps.value();
}

//void class_state_infinite::assert_positions() const {
//    if(not math::all_equal(Lblock->get_position(), Lblock2->get_position(), HA->get_position(), MPS->MPS_A->get_position()))
//        throw std::runtime_error(fmt::format("Position error: Lblock ({}), Lblock2 ({}), HA ({}), MPS_A ({})", Lblock->get_position(), Lblock2->get_position(),
//                                             HA->get_position(), MPS->MPS_A->get_position()));
//
//    if(not math::all_equal(Rblock->get_position(), Rblock2->get_position(), HB->get_position(), MPS->MPS_B->get_position()))
//        throw std::runtime_error(fmt::format("Position error: Rblock ({}), Rblock2 ({}), HB ({}), MPS_B ({})", Rblock->get_position(), Rblock2->get_position(),
//                                             HB->get_position(), MPS->MPS_B->get_position()));
//}

//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//

// Eigen::Tensor<Scalar,4> class_state_infinite::optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, eigutils::eigSetting::Ritz ritz){
//    std::array<long,4> shape_theta4 = theta.dimensions();
//    std::array<long,4> shape_mpo4   = HA->MPO().dimensions();
//
//    t_eig.tic();
//    int nev = 1;
//    using namespace settings::precision;
//    using namespace eigutils::eigSetting;
//    DenseHamiltonianProduct<Scalar>  matrix (Lblock->block.data(), Rblock->block.data(), HA->MPO().data(), HB->MPO().data(), shape_theta4, shape_mpo4);
//    class_eigsolver solver;
//    solver.eigs_dense(matrix,nev,eig_max_ncv,NAN,Form::SYMMETRIC,ritz,Side::R,true,true);
//
//    auto eigvals           = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data()
//    ,solver.solution.meta.cols); auto eigvecs           = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX,
//    Form::SYMMETRIC>().data(),solver.solution.meta.rows);
//
//    t_eig.toc();
//    t_eig.print_last_time_interval();
//
//    E_optimal = std::real(eigvals(0));
//    return eigvecs.reshape(theta.dimensions());
//}

//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//
// Eigen::Tensor<Scalar, 4> class_state_infinite::evolve_MPS(const Eigen::Tensor<Scalar, 4> &U)
///*!
//@verbatim
//  1--[ Θ ]--3
//     |   |
//     0   2
//                   1--[ Θ ]--3
//     0   1   --->     |   |
//     |   |            0   2
//     [ U ]
//     |   |
//     2   3
//@endverbatim
//*/
//
//{
//    return U.contract(MPS->get_theta(), idx({0,1},{0,2}))
//            .shuffle(array4{0,2,1,3});
//}
//
// Eigen::Tensor<Scalar, 4> class_state_infinite::evolve_MPS(const Eigen::Tensor<Scalar, 4> &theta, const Eigen::Tensor<Scalar, 4> &U)
///*!
//@verbatim
//  1--[ Θ ]--3
//     |   |
//     0   2
//                   1--[ Θ ]--3
//     0   1   --->     |   |
//     |   |            0   2
//     [ U ]
//     |   |
//     2   3
//@endverbatim
//*/
//{
//    return U.contract(theta, idx({0,1},{0,2}))
//            .shuffle(array4{0,2,1,3});
//}

//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS->
//============================================================================//
// Eigen::Tensor<Scalar,4> class_state_infinite::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta,long chi_, double svd_threshold){
//    class_SVD SVD;
//    SVD.setThreshold(svd_threshold);
//    auto[U, S, V] = SVD.schmidt(theta,chi_);
//    MPS->truncation_error         = SVD.get_truncation_error();
//    MPS->LC  = S;
//    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
//    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS->MPS_B->get_L()), idx({2},{0}));
//    MPS->MPS_A->set_G(L_U);
//    MPS->MPS_B->set_G(V_L);
//    return get_theta();
//}
//
// void class_state_infinite::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta, class_mps_2site &MPS_out,long chi_, double svd_threshold){
//    class_SVD SVD;
//    SVD.setThreshold(svd_threshold);
//    auto[U, S, V] = SVD.schmidt(theta, chi_);
//    MPS_out.truncation_error = SVD.get_truncation_error();
//    MPS_out.LC  = S;
//    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS_out.MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
//    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS_out.MPS_B->get_L()), idx({2},{0}));
//    MPS_out.MPS_A->set_G(L_U);
//    MPS_out.MPS_B->set_G(V_L);
//}
//

void class_state_infinite::assert_validity() const { MPS->assert_validity(); }
bool class_state_infinite::is_real() const { return MPS->is_real(); }
bool class_state_infinite::has_nan() const { return MPS->has_nan(); }
//
//template<typename T>
//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_matrix() const {
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
//    long                shape   = MPS->chiA() * MPS->spindim() * MPS->chiB() * MPS->spin_dim_A();
//    Eigen::Tensor<T, 8> H_local = tempL.contract(tempR, Textra::idx({4}, {4})).shuffle(Textra::array8{0, 1, 4, 5, 2, 3, 6, 7});
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(H_local.data(), shape, shape);
//}

//template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>               class_state_infinite::get_H_local_matrix<double>() const;
//template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_matrix<std::complex<double>>() const;
//
//template<typename T>
//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_sq_matrix() const {
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
//    long                shape   = MPS->chiA() * MPS->spindim() * MPS->chiB() * MPS->spin_dim_A();
//    Eigen::Tensor<T, 8> H_local = tempL.contract(tempR, Textra::idx({4, 5}, {4, 5})).shuffle(Textra::array8{0, 1, 4, 5, 2, 3, 6, 7});
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(H_local.data(), shape, shape);
//}

//template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>               class_state_infinite::get_H_local_sq_matrix<double>() const;
//template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> class_state_infinite::get_H_local_sq_matrix<std::complex<double>>() const;

// void class_state_infinite::enlarge_environment(int direction) {
//    assert_positions();
//    auto position_A_new = MPS->MPS_A->get_position() + 1;
//    *Lblock             = Lblock->enlarge(*MPS->MPS_A, *HA);
//    *Rblock             = Rblock->enlarge(*MPS->MPS_B, *HB);
//    *Lblock2            = Lblock2->enlarge(*MPS->MPS_A, *HA);
//    *Rblock2            = Rblock2->enlarge(*MPS->MPS_B, *HB);
//    set_positions(position_A_new);
//    if(direction != 0) throw std::runtime_error("Ooops, direction != 0");
//    if (direction == 1){
//        *Lblock  = Lblock->enlarge(*MPS->MPS_A,  *HA);
//        *Lblock2 = Lblock2->enlarge(*MPS->MPS_A, *HA);
//        Lblock->set_position (HB->get_position());
//        Lblock2->set_position(HB->get_position());
//    }else if (direction == -1){
//        *Rblock  = Rblock->enlarge(*MPS->MPS_B,  *HB);
//        *Rblock2 = Rblock2->enlarge(*MPS->MPS_B, *HB);
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
    MPS->MPS_A->set_position(position);
    MPS->MPS_B->set_position(position + 1);
}

void class_state_infinite::set_mps(const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB) {
    MPS->set_mps(MA, LC, MB);
}
void class_state_infinite::set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC,
                                   const Eigen::Tensor<Scalar, 3> &MB, const Eigen::Tensor<Scalar, 1> &LB) {
    MPS->set_mps(LA, MA, LC, MB, LB);
}

void class_state_infinite::clear_cache() const { cache = Cache(); }

void class_state_infinite::clear_measurements() const {
    measurements                              = state_measure_infinite();
    tools::common::views::components_computed = false;
}

void class_state_infinite::do_all_measurements() const { tools::infinite::measure::do_all_measurements(*this); }

void class_state_infinite::swap_AB() { MPS->swap_AB(); }
