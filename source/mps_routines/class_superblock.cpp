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
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_superblock::~class_superblock()=default;



size_t class_superblock::get_length() const {
    size_t length = Lblock->size + Rblock->size + 2;
    assert(length == environment_size + 2 and "ERROR: System length mismatch!");
    return length;

}


//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//



Eigen::Tensor<Scalar,4> class_superblock::optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, eigsolver_properties::Ritz ritz){
    std::array<long,4> shape_theta4 = theta.dimensions();
    std::array<long,4> shape_mpo4   = HA->MPO.dimensions();

    t_eig.tic();
    int nev = 1;
    using namespace settings::precision;
    class_eigsolver_arpack<Scalar, eigsolver_properties::Form::GENERAL> solver(eigThreshold,eigMaxIter,eigMaxNcv, true, false);
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



Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_matrix (){
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,5>tempL   = Lblock->block.contract(HA->MPO,Textra::idx({2},{0})).shuffle(Textra::array5{4,1,3,0,2});
    Eigen::Tensor<Scalar,5>tempR   = Rblock->block.contract(HB->MPO,Textra::idx({2},{1})).shuffle(Textra::array5{4,1,3,0,2});
    Eigen::Tensor<Scalar,8>H_local = tempL.contract(tempR, Textra::idx({4},{4})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//    tempL.resize(0,0,0,0,0);
//    tempR.resize(0,0,0,0,0);
//    tempH = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7}).eval();
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_matrix_real (){
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
//    std::cout << "STATUS: Computing tempL \n";
    Eigen::Tensor<double,5>tempL   = Lblock->block.contract(HA->MPO,Textra::idx({2},{0})).real().shuffle(Textra::array5{4,1,3,0,2});
//    std::cout << "STATUS: tempL dims: " << tempL.dimensions() << std::endl;
//    std::cout << "STATUS: Computing tempR \n";
    Eigen::Tensor<double,5>tempR   = Rblock->block.contract(HB->MPO,Textra::idx({2},{1})).real().shuffle(Textra::array5{4,1,3,0,2});
//    std::cout << "STATUS: tempR dims: " << tempR.dimensions() << std::endl;
//    std::cout << "STATUS: Computing H_local \n";
    Eigen::Tensor<double,8>H_local = tempL.contract(tempR, Textra::idx({4},{4})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//    tempL.resize(0,0,0,0,0);
//    tempR.resize(0,0,0,0,0);
//    tempH = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7}).eval();
//    Eigen::Tensor<double,8> H_local = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7}).eval();
//    std::cout << "STATUS: Returning get_H_local_matrix_real\n";
    return Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_sq_matrix (){
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
                                                 .contract(HA->MPO,Textra::idx({2,5},{0,2}))
                                                 .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
                                                 .contract(HB->MPO,Textra::idx({2,5},{1,2}))
                                                 .shuffle(Textra::array6{5,1,3,0,2,4});

    Eigen::Tensor<Scalar,8> H_local = tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//    tempL.resize(0,0,0,0,0,0);
//    tempR.resize(0,0,0,0,0,0);
//    Eigen::Tensor<Scalar,8> H_local = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//    tempH = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7}).eval();
//
//
//    Eigen::Tensor<Scalar,8> tempH =  Lblock2->block
//            .contract(HA->MPO      , Textra::idx({2},{0}))
//            .contract(HA->MPO      , Textra::idx({2,5},{0,2}))
//            .contract(HB->MPO      , Textra::idx({2},{0}))
//            .contract(HB->MPO      , Textra::idx({3,7},{0,2}))
//            .contract(Rblock2->block, Textra::idx({4,6},{2,3}))
//            .shuffle(Textra::array8{3,1,5,7,2,0,4,6});
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}

Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_sq_matrix_real (){
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<double,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
                                                 .contract(HA->MPO,Textra::idx({2,5},{0,2}))
                                                 .real().shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<double,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
                                                 .contract(HB->MPO,Textra::idx({2,5},{1,2}))
                                                 .real().shuffle(Textra::array6{5,1,3,0,2,4});

    Eigen::Tensor<double,8>H_local = tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//    tempL.resize(0,0,0,0,0,0);
//    tempR.resize(0,0,0,0,0,0);
//    tempH = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7}).eval();
    return Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


Textra::SparseMatrixType<Scalar> class_superblock::get_H_local_sparse_matrix (double prune){
    return Textra::Tensor2_to_SparseMatrix(get_H_local_rank2(),prune);
}
Textra::SparseMatrixType<Scalar> class_superblock::get_H_local_sq_sparse_matrix (double prune){
    return Textra::Tensor2_to_SparseMatrix(get_H_local_sq_rank2(),prune);
}


Eigen::Tensor<Scalar,2> class_superblock::get_H_local_rank2 (){
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7}).reshape(Textra::array2{shape, shape});

//
//    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
//    return Lblock->block
//            .contract(HA->MPO      , Textra::idx({2},{0}))
//            .contract(HB->MPO      , Textra::idx({2},{0}))
//            .contract(Rblock->block, Textra::idx({4},{2}))
//            .shuffle(Textra::array8{3,1,5,7,2,0,4,6})
//            .reshape(Textra::array2{shape, shape});
}

Eigen::Tensor<Scalar,8> class_superblock::get_H_local_rank8 (){
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//
//
//    return Lblock->block
//            .contract(HA->MPO      , Textra::idx({2},{0}))
//            .contract(HB->MPO      , Textra::idx({2},{0}))
//            .contract(Rblock->block, Textra::idx({4},{2}))
//            .shuffle(Textra::array8{3,1,5,7,2,0,4,6});
}

Eigen::Tensor<Scalar,2> class_superblock::get_H_local_sq_rank2 (){

    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7}).reshape(Textra::array2{shape, shape});

//
//    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
//    return Lblock2->block
//            .contract(HA->MPO      , Textra::idx({2},{0}))
//            .contract(HA->MPO      , Textra::idx({2,5},{0,2}))
//            .contract(HB->MPO      , Textra::idx({2},{0}))
//            .contract(HB->MPO      , Textra::idx({3,7},{0,2}))
//            .contract(Rblock2->block, Textra::idx({4,6},{2,3}))
//            .shuffle(Textra::array8{3,1,5,7,2,0,4,6})
//            .reshape(Textra::array2{shape, shape});
}

Eigen::Tensor<Scalar,8> class_superblock::get_H_local_sq_rank8 (){
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//
//    return Lblock2->block
//            .contract(HA->MPO      , Textra::idx({2},{0}))
//            .contract(HA->MPO      , Textra::idx({2,5},{0,2}))
//            .contract(HB->MPO      , Textra::idx({2},{0}))
//            .contract(HB->MPO      , Textra::idx({3,7},{0,2}))
//            .contract(Rblock2->block, Textra::idx({4,6},{2,3}))
//            .shuffle(Textra::array8{3,1,5,7,2,0,4,6});
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
//        assert(Lblock->get_position()  == HA->get_position());
//        assert(Lblock2->get_position() == HA->get_position());
//        assert(Rblock->get_position()  == HB->get_position());
//        assert(Rblock2->get_position() == HB->get_position());
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


void class_superblock::set_superblock(
        const Eigen::Tensor<Scalar,4> & Lblock2_,
        const Eigen::Tensor<Scalar,3> & Lblock_,
        const Eigen::Tensor<Scalar,4> & MPO_A,
        const Eigen::Tensor<Scalar,1> & LA,
        const Eigen::Tensor<Scalar,3> & GA,
        const Eigen::Tensor<Scalar,1> & LC,
        const Eigen::Tensor<Scalar,3> & GB,
        const Eigen::Tensor<Scalar,1> & LB,
        const Eigen::Tensor<Scalar,4> & MPO_B,
        const Eigen::Tensor<Scalar,3> & Rblock_,
        const Eigen::Tensor<Scalar,4> & Rblock2_
        )
{
    Lblock2->block  = Lblock2_;
    Lblock->block   = Lblock_;
    HA->MPO = MPO_A;
    MPS->MPS_A->set_L(LA);
    MPS->MPS_A->set_G(GA);
    MPS->LC = LC;
    MPS->MPS_B->set_G(GB);
    MPS->MPS_B->set_L(LB);
    HB->MPO = MPO_B;
    Rblock ->block = Rblock_;
    Rblock2->block = Rblock2_;
}



void  class_superblock::swap_AB(){
    MPS->swap_AB();
}

