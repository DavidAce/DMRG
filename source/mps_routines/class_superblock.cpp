//
// Created by david on 7/22/17.
//

//#include <mps_routines/class_optimize_mps.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_hamiltonian.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <iomanip>
//#include <SymEigsSolver.h>
//#include <MatOp/SparseSymMatProd.h>
#include <general/class_svd_wrapper.h>
//#include <mps_routines/class_custom_contraction.h>
#include <sim_parameters/nmspc_sim_settings.h>
//#include <arpackpp/arscomp.h>
#include <general/class_arpack_eigsolver.h>

#define profile_optimization 0

using namespace std;
using namespace Textra;
using Scalar = class_superblock::Scalar;

class_superblock::class_superblock():
        MPS(std::make_shared<class_mps>()),
        H(std::make_shared<class_mpo>()),
        HA(std::make_shared<class_hamiltonian>()),
        HB(std::make_shared<class_hamiltonian>()),
        Lblock(std::make_shared<class_environment>("L")),
        Rblock(std::make_shared<class_environment>("R")),
        Lblock2(std::make_shared<class_environment_var>("L")),
        Rblock2(std::make_shared<class_environment_var>("R")),
        SVD(std::make_shared<class_SVD<Scalar>>())
{
    t_eig.set_properties(profile_optimization, 10,"Time optimizing ");
    MPS->initialize(H->local_dimension);
    HA->set_parameters(settings::model::J, settings::model::g, 0.0);
    HB->set_parameters(settings::model::J, settings::model::g, 0.0);
    set_current_dimensions();

}


//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//



Textra::Tensor<Scalar,4> class_superblock::optimize_MPS(Textra::Tensor<Scalar, 4> &theta){
    std::array<long,4> shape_theta4 = theta.dimensions();
    std::array<long,4> shape_mpo4   = HA->MPO.dimensions();

    t_eig.tic();
    int nev = 1;
    using namespace settings::precision;
    class_arpack_eigsolver<Scalar, Form::GENERAL> solver(Lblock->block.data(), Rblock->block.data(), HA->MPO.data(), HB->MPO.data(), shape_theta4, shape_mpo4, Ritz::SR, nev, eigMaxNcv,eigThreshold,eigMaxIter, true ,theta.data());

    TensorMap<const Tensor<const Scalar,2>> eigvecs (solver.ref_eigvecs().data(), shape_theta4[0]*shape_theta4[1], shape_theta4[2]*shape_theta4[3]);
    TensorMap<const Tensor<const Scalar,1>> eigvals (solver.ref_eigvals().data(), nev);
    t_eig.toc();
    t_eig.print_delta();

    E_one_site = std::real(eigvals(0))/2.0;
//    double L = Lblock->size + Rblock->size;
//    std::cout <<setprecision(16) << "E_lanczos: " <<  std::real(eigvals(0))  << " L : " << L << " " << environment_size + 2 << std::endl;
    //    using namespace chrono;
//    std::cout << "Time: " << duration_cast<duration<double>>(t_eig.measured_time).count() << " ";
//    t_eig.print_delta();
//    std::cout << " iter: " << opt.iter << " counter: " << opt.counter << " E_one_site: " << E_one_site <<  " shape4: " << theta.dimensions() << "\n";
    return eigvecs.reshape(shape_theta4);
}




//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//
Textra::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Textra::Tensor<Scalar, 4> &U)
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

Textra::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Textra::Tensor<Scalar, 4> &theta, const Textra::Tensor<Scalar, 4> &U)
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
Textra::Tensor<Scalar,4> class_superblock::truncate_MPS(const Textra::Tensor<Scalar, 4> &theta,long chi_max_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta, d, chiL, chi_max_, chiR);
    MPS->truncation_error = SVD->get_truncation_error();
    MPS->LA  = S;
    MPS->GA  = asDiagonalInversed(MPS->LB_left).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    MPS->GB  = V.contract(asDiagonalInversed(MPS->LB), idx({2},{0}));
    return MPS->get_theta();
}

void class_superblock::truncate_MPS(const Textra::Tensor<Scalar, 4> &theta, const std::shared_ptr<class_mps> &MPS_out,long chi_max_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta, d, chiL, chi_max_, chiR);
    MPS_out->truncation_error = SVD->get_truncation_error();
    MPS_out->LA  = S;
    MPS_out->GA  = asDiagonalInversed(MPS_out->LB_left).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    MPS_out->GB  = V.contract(asDiagonalInversed(MPS_out->LB), idx({2},{0}));
}



void class_superblock::enlarge_environment(int direction){
    if (direction == 1){
        Lblock->enlarge(MPS,  HA->MPO_reduced());
        Lblock2->enlarge(MPS, HA->MPO_reduced());
    }else if (direction == -1){
        Rblock->enlarge(MPS,  HB->MPO_reduced());
        Rblock2->enlarge(MPS, HB->MPO_reduced());
    }else if(direction == 0){
        Lblock->enlarge(MPS,  HA->MPO_reduced());
        Rblock->enlarge(MPS,  HB->MPO_reduced());
        Lblock2->enlarge(MPS, HA->MPO_reduced());
        Rblock2->enlarge(MPS, HB->MPO_reduced());
        environment_size = Lblock->size + Rblock->size;
    }
}


void class_superblock::set_current_dimensions(){
    d       = H->local_dimension;
    chi     = MPS->GA.dimension(2);
    chiL    = MPS->GA.dimension(1);
    chiR    = MPS->GB.dimension(2);
    shape1  = {d * chiL * d * chiR};
    shape2  = {d * chiL , d * chiR};
    shape4  = {d , chiL , d , chiR};
}


void  class_superblock::swap_AB(){
    MPS->swap_AB();
}

