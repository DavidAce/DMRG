//
// Created by david on 7/22/17.
//

//#include <mps_routines/class_optimize_mps.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps_2site.h>
#include <iomanip>
#include <general/class_svd_wrapper.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/class_eigsolver_arpack.h>
#include <model/class_hamiltonian_factory.h>
#define profile_optimization 0

using namespace std;
using namespace Textra;
using Scalar = class_superblock::Scalar;

class_superblock::class_superblock():
        MPS(std::make_unique<class_mps_2site>()),
        HA (class_hamiltonian_factory::create_mpo(settings::model::model_type)),
        HB (class_hamiltonian_factory::create_mpo(settings::model::model_type)),
        Lblock(std::make_unique<class_environment>("L")),
        Rblock(std::make_unique<class_environment>("R")),
        Lblock2(std::make_unique<class_environment_var>("L")),
        Rblock2(std::make_unique<class_environment_var>("R")),
        SVD(std::make_unique<class_SVD<Scalar>>())
{
    HA->set_position(0);
    HB->set_position(1);
    t_eig.set_properties(profile_optimization, 10,"Time optimizing ");
    spin_dimension = HA->get_spin_dimension();
    MPS->initialize(spin_dimension);
}




// We need to make a destructor manually for the enclosing class "class_superblock"
// that encloses "class_hamiltonian_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_hamiltonian_base"
// Read mode: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_superblock::~class_superblock()=default;

//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//



Eigen::Tensor<Scalar,4> class_superblock::optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, Ritz ritz){
    std::array<long,4> shape_theta4 = theta.dimensions();
    std::array<long,4> shape_mpo4   = HA->MPO.dimensions();

    t_eig.tic();
    int nev = 1;
    using namespace settings::precision;
    class_eigsolver_arpack<Scalar, Form::GENERAL> solver(eigThreshold,eigMaxIter,eigMaxNcv, true, false);
    solver.optimize_mps(Lblock->block.data(), Rblock->block.data(), HA->MPO.data(), HB->MPO.data(), shape_theta4, shape_mpo4, nev, eigMaxNcv, ritz , false ,theta.data());

    Eigen::TensorMap<const Eigen::Tensor<const Scalar,2>> eigvecs (solver.ref_eigvecs().data(), shape_theta4[0]*shape_theta4[1], shape_theta4[2]*shape_theta4[3]);
    Eigen::TensorMap<const Eigen::Tensor<const Scalar,1>> eigvals (solver.ref_eigvals().data(), nev);
    t_eig.toc();
    t_eig.print_delta();

    E_optimal = std::real(eigvals(0));
//    double L = Lblock->size + Rblock->size;
//    std::cout <<setprecision(16) << "E_lanczos: " <<  std::real(eigvals(0))  << " L : " << L << " " << environment_size + 2 << std::endl;
    //    using namespace chrono;
//    std::cout << "Time: " << duration_cast<duration<double>>(t_eig.measured_time).count() << " ";
//    t_eig.print_delta();
//    std::cout << " iter: " << opt.iter << " counter: " << opt.counter << " E_optimal: " << E_optimal <<  " shape4: " << theta.dimensions() << "\n";
    return eigvecs.reshape(shape_theta4);
}




//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//
Eigen::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Eigen::Tensor<Scalar, 4> &U)
/*!
@verbatim
  1--[ Θ ]--3
     |   |
     0   2
                   1--[ Θ ]--3
     0   1   --->     |   |
     |   |            0   2
     [ U ]
     |   |
     2   3
@endverbatim
*/

{
    return U.contract(MPS->get_theta(), idx({0,1},{0,2}))
            .shuffle(array4{0,2,1,3});
}

Eigen::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Eigen::Tensor<Scalar, 4> &theta, const Eigen::Tensor<Scalar, 4> &U)
/*!
@verbatim
  1--[ Θ ]--3
     |   |
     0   2
                   1--[ Θ ]--3
     0   1   --->     |   |
     |   |            0   2
     [ U ]
     |   |
     2   3
@endverbatim
*/
{
    return U.contract(theta, idx({0,1},{0,2}))
            .shuffle(array4{0,2,1,3});
}

//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS->
//============================================================================//
Eigen::Tensor<Scalar,4> class_superblock::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta,long chi_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta,chi_);
    MPS->truncation_error = SVD->get_truncation_error();
    MPS->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS->MPS_B->get_L()), idx({2},{0}));
    MPS->MPS_A->set_G(L_U);
    MPS->MPS_B->set_G(V_L);
//    std::cout << setprecision(3) << "S: \n "<< S << "\nU\n" << U << "\nV\n" << V<< std::endl;
    return MPS->get_theta();
}

void class_superblock::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta, const std::shared_ptr<class_mps_2site> &MPS_out,long chi_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta, chi_);
    MPS_out->truncation_error = SVD->get_truncation_error();
    MPS_out->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS_out->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS_out->MPS_B->get_L()), idx({2},{0}));
    MPS_out->MPS_A->set_G(L_U);
    MPS_out->MPS_B->set_G(V_L);

}



void class_superblock::enlarge_environment(int direction){
    if (direction == 1){
        assert(Lblock->get_position()  == HA->get_position());
        assert(Lblock2->get_position() == HA->get_position());
        Lblock->enlarge(MPS,  HA->MPO_reduced_view());
        Lblock2->enlarge(MPS, HA->MPO_reduced_view());
        Lblock->set_position (HB->get_position());
        Lblock2->set_position(HB->get_position());
    }else if (direction == -1){
        assert(Rblock->get_position()  == HB->get_position());
        assert(Rblock2->get_position() == HB->get_position());
        Rblock->enlarge(MPS,  HB->MPO_reduced_view());
        Rblock2->enlarge(MPS, HB->MPO_reduced_view());
        Rblock->set_position (HA->get_position());
        Rblock2->set_position(HA->get_position());
    }else if(direction == 0){
        assert(Lblock->get_position()  == HA->get_position());
        assert(Lblock2->get_position() == HA->get_position());
        assert(Rblock->get_position()  == HB->get_position());
        assert(Rblock2->get_position() == HB->get_position());
        Lblock->enlarge(MPS,  HA->MPO_reduced_view());
        Rblock->enlarge(MPS,  HB->MPO_reduced_view());
        Lblock2->enlarge(MPS, HA->MPO_reduced_view());
        Rblock2->enlarge(MPS, HB->MPO_reduced_view());
        Lblock->set_position (HB->get_position());
        Rblock->set_position (HB->get_position()+1);
        Lblock2->set_position(HB->get_position());
        Rblock2->set_position(HB->get_position()+1);
        environment_size = Lblock->size + Rblock->size;
    }
}



void  class_superblock::swap_AB(){
    MPS->swap_AB();
}

