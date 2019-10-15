//
// Created by david on 7/22/17.
//

//#include <state/class_optimize_mps.h>
#include <state/class_infinite_state.h>
#include <state/class_environment.h>
#include <state/class_mps_site.h>
#include <state/class_mps_2site.h>
#include <tools/nmspc_tools.h>
#include <iomanip>
#include <math/class_eigsolver.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <simulation/nmspc_settings.h>
#include <model/class_model_factory.h>
#define profile_optimization 0

using namespace std;
using namespace Textra;
using Scalar = class_infinite_state::Scalar;

class_infinite_state::class_infinite_state(SimulationType sim_type_,std::string sim_name_):
        sim_type(sim_type_),
        sim_name(sim_name_),
        MPS(std::make_unique<class_mps_2site>()),
        HA (class_model_factory::create_mpo(0,settings::model::model_type)),
        HB (class_model_factory::create_mpo(1,settings::model::model_type)),
        Lblock(std::make_unique<class_environment>("L")),
        Rblock(std::make_unique<class_environment>("R")),
        Lblock2(std::make_unique<class_environment_var>("L")),
        Rblock2(std::make_unique<class_environment_var>("R"))
{
    log = Logger::setLogger(sim_name,settings::console::verbosity);
//    t_eig.set_properties(profile_optimization, 10,"Time optimizing ");
    MPS->initialize(HA->get_spin_dimension());
}






// We need to make a destructor manually for the enclosing class "class_infinite_state"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
//class_infinite_state::~class_infinite_state()=default;

// we could the default move constructor
//class_infinite_state::class_infinite_state(class_infinite_state&&)  = default;
// we could use the default move assignment operator
//class_infinite_state & class_infinite_state::operator=(class_infinite_state&&) = default;

// ...but we would have to provide a deep copying assignment operator
//class_infinite_state & class_infinite_state::operator=(class_infinite_state const& source) {
//    HA = source.HA->clone();
//    HB = source.HB->clone();
//    return *this;
//}

void class_infinite_state::clear(){
    log->trace("Copying state to state");
    *this = class_infinite_state(sim_type,sim_name);
}

size_t class_infinite_state::get_length() const {
    size_t length = Lblock->sites + Rblock->sites + 2;
    return length;
}

size_t class_infinite_state::get_position() const { return Lblock->sites;}

long class_infinite_state::get_chi() const {return MPS->chiC();}

long class_infinite_state::get_chi_lim()  const {
    //Should get the the current limit on allowed bond dimension
    return chi_lim;
}
void class_infinite_state::set_chi_lim(long chi_lim_){
    //Should set the the current limit on allowed bond dimension
    chi_lim = chi_lim_;
}

Eigen::DSizes<long,4> class_infinite_state::dimensions()const{return MPS->dimensions();}

Eigen::Tensor<Scalar, 4> class_infinite_state::get_theta() const { return MPS->get_theta();}



//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//



//Eigen::Tensor<Scalar,4> class_infinite_state::optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, eigutils::eigSetting::Ritz ritz){
//    std::array<long,4> shape_theta4 = theta.dimensions();
//    std::array<long,4> shape_mpo4   = HA->MPO().dimensions();
//
//    t_eig.tic();
//    int nev = 1;
//    using namespace settings::precision;
//    using namespace eigutils::eigSetting;
//    DenseHamiltonianProduct<Scalar>  matrix (Lblock->block.data(), Rblock->block.data(), HA->MPO().data(), HB->MPO().data(), shape_theta4, shape_mpo4);
//    class_eigsolver solver;
//    solver.eigs_dense(matrix,nev,eigMaxNcv,NAN,Form::SYMMETRIC,ritz,Side::R,true,true);
//
//    auto eigvals           = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
//    auto eigvecs           = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);
//
//    t_eig.toc();
//    t_eig.print_delta();
//
//    E_optimal = std::real(eigvals(0));
//    return eigvecs.reshape(theta.dimensions());
//}




//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//
//Eigen::Tensor<Scalar, 4> class_infinite_state::evolve_MPS(const Eigen::Tensor<Scalar, 4> &U)
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
//Eigen::Tensor<Scalar, 4> class_infinite_state::evolve_MPS(const Eigen::Tensor<Scalar, 4> &theta, const Eigen::Tensor<Scalar, 4> &U)
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
//Eigen::Tensor<Scalar,4> class_infinite_state::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta,long chi_, double SVDThreshold){
//    class_SVD SVD;
//    SVD.setThreshold(SVDThreshold);
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
//void class_infinite_state::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta, class_mps_2site &MPS_out,long chi_, double SVDThreshold){
//    class_SVD SVD;
//    SVD.setThreshold(SVDThreshold);
//    auto[U, S, V] = SVD.schmidt(theta, chi_);
//    MPS_out.truncation_error = SVD.get_truncation_error();
//    MPS_out.LC  = S;
//    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS_out.MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
//    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS_out.MPS_B->get_L()), idx({2},{0}));
//    MPS_out.MPS_A->set_G(L_U);
//    MPS_out.MPS_B->set_G(V_L);
//}
//



bool class_infinite_state::isReal() const {
    bool MPS_isReal         = MPS->isReal();
    bool HA_isReal          = HA->isReal();
    bool HB_isReal          = HB->isReal();
    bool Lblock_isReal      = Lblock->isReal();
    bool Rblock_isReal      = Rblock->isReal();
    bool Lblock2_isReal     = Lblock2->isReal();
    bool Rblock2_isReal     = Rblock2->isReal();

    return
        MPS_isReal
    and HA_isReal
    and HB_isReal
    and Lblock_isReal
    and Rblock_isReal
    and Lblock2_isReal
    and Rblock2_isReal;
}

template<typename T>
Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> class_infinite_state::get_H_local_matrix ()const{
    Eigen::Tensor<T,5>tempL;
    Eigen::Tensor<T,5>tempR;
    if constexpr (std::is_same<T,double>::value){
        if(not Lblock->isReal()){throw std::runtime_error("Discarding imaginary data from Lblock when building H_local");}
        if(not Rblock->isReal()){throw std::runtime_error("Discarding imaginary data from Rblock when building H_local");}
        if(not HA->isReal())    {throw std::runtime_error("Discarding imaginary data from MPO A when building H_local");}
        if(not HB->isReal())    {throw std::runtime_error("Discarding imaginary data from MPO B when building H_local");}
        tempL   = Lblock->block.contract(HA->MPO(),Textra::idx({2},{0})).real().shuffle(Textra::array5{4,1,3,0,2}).real();
        tempR   = Rblock->block.contract(HB->MPO(),Textra::idx({2},{1})).real().shuffle(Textra::array5{4,1,3,0,2}).real();
    }else{
        tempL   = Lblock->block.contract(HA->MPO(),Textra::idx({2},{0})).shuffle(Textra::array5{4,1,3,0,2});
        tempR   = Rblock->block.contract(HB->MPO(),Textra::idx({2},{1})).shuffle(Textra::array5{4,1,3,0,2});
    }
    long shape = MPS->chiA() * MPS->spindim() * MPS->chiB() * MPS->spindim();
    Eigen::Tensor<T,8>H_local = tempL.contract(tempR, Textra::idx({4},{4})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}

template Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>                class_infinite_state::get_H_local_matrix<double>() const;
template Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic>  class_infinite_state::get_H_local_matrix<std::complex<double>>() const;


template<typename T>
Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> class_infinite_state::get_H_local_sq_matrix ()const{
    Eigen::Tensor<T,6>tempL;
    Eigen::Tensor<T,6>tempR;
    if constexpr(std::is_same<T,double>::value){
        if(not Lblock2->isReal()){throw std::runtime_error("Discarding imaginary data from Lblock2 when building H_local_sq");}
        if(not Rblock2->isReal()){throw std::runtime_error("Discarding imaginary data from Rblock2 when building H_local_sq");}
        if(not HA->isReal())     {throw std::runtime_error("Discarding imaginary data from MPO A when building H_local_sq");}
        if(not HB->isReal())     {throw std::runtime_error("Discarding imaginary data from MPO B when building H_local_sq");}
        tempL = Lblock2->block
                .contract(HA->MPO(),Textra::idx({2},{0}))
                .contract(HA->MPO(),Textra::idx({2,5},{0,2}))
                .real()
                .shuffle(Textra::array6{5,1,3,0,2,4});
        tempR = Rblock2->block
                .contract(HB->MPO(),Textra::idx({2},{1}))
                .contract(HB->MPO(),Textra::idx({2,5},{1,2}))
                .real()
                .shuffle(Textra::array6{5,1,3,0,2,4});
    }else{
        tempL = Lblock2->block
                .contract(HA->MPO(),Textra::idx({2},{0}))
                .contract(HA->MPO(),Textra::idx({2,5},{0,2}))
                .shuffle(Textra::array6{5,1,3,0,2,4});
        tempR = Rblock2->block
                .contract(HB->MPO(),Textra::idx({2},{1}))
                .contract(HB->MPO(),Textra::idx({2,5},{1,2}))
                .shuffle(Textra::array6{5,1,3,0,2,4});
    }

    long shape = MPS->chiA() * MPS->spindim() * MPS->chiB() * MPS->spindim();
    Eigen::Tensor<T,8> H_local = tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


template Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>                class_infinite_state::get_H_local_sq_matrix<double>() const ;
template Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic>  class_infinite_state::get_H_local_sq_matrix<std::complex<double>>() const;



void class_infinite_state::enlarge_environment(int direction){
    if (direction == 1){
        assert(Lblock->get_position()  == HA->get_position());
        assert(Lblock2->get_position() == HA->get_position());
        Lblock->enlarge(*MPS->MPS_A,  HA->MPO_reduced_view());
        Lblock2->enlarge(*MPS->MPS_A, HA->MPO_reduced_view());
        Lblock->set_position (HB->get_position());
        Lblock2->set_position(HB->get_position());
    }else if (direction == -1){
        assert(Rblock->get_position()  == HB->get_position());
        assert(Rblock2->get_position() == HB->get_position());
        Rblock->enlarge(*MPS->MPS_B,  HB->MPO_reduced_view());
        Rblock2->enlarge(*MPS->MPS_B, HB->MPO_reduced_view());
        Rblock->set_position (HA->get_position());
        Rblock2->set_position(HA->get_position());
    }else if(direction == 0){
        Lblock->enlarge (*MPS->MPS_A,  HA->MPO_reduced_view());
        Rblock->enlarge (*MPS->MPS_B,  HB->MPO_reduced_view());
        Lblock2->enlarge(*MPS->MPS_A, HA->MPO_reduced_view());
        Rblock2->enlarge(*MPS->MPS_B, HB->MPO_reduced_view());

        Lblock->set_position (HB->get_position());
        Rblock->set_position (HB->get_position()+1);
        Lblock2->set_position(HB->get_position());
        Rblock2->set_position(HB->get_position()+1);
    }
    unset_measurements();
}


void class_infinite_state::set_positions(int position){
    MPS->MPS_A->set_position(position);
    MPS->MPS_B->set_position(position+1);
    Lblock->set_position(position);
    Rblock->set_position(position+1);
    Lblock2->set_position(position);
    Rblock2->set_position(position+1);
    HA->set_position(position);
    HB->set_position(position+1);
}


void class_infinite_state::unset_measurements() const {
    measurements = Measurements();
    tools::common::views::components_computed = false;
}

void class_infinite_state::do_all_measurements() const {
    measurements.length                         = tools::infinite::measure::length(*this);
    measurements.bond_dimension                 = tools::infinite::measure::bond_dimension(*this);
    measurements.norm                           = tools::infinite::measure::norm(*this);
    measurements.truncation_error               = tools::infinite::measure::truncation_error(*this);
    measurements.energy_mpo                     = tools::infinite::measure::energy_mpo(*this);  //This number is needed for variance calculation!
    measurements.energy_per_site_mpo            = tools::infinite::measure::energy_per_site_mpo(*this);
    measurements.energy_per_site_ham            = tools::infinite::measure::energy_per_site_ham(*this);
    measurements.energy_per_site_mom            = tools::infinite::measure::energy_per_site_mom(*this);
    measurements.energy_variance_mpo            = tools::infinite::measure::energy_variance_mpo(*this);
    measurements.energy_variance_per_site_mpo   = tools::infinite::measure::energy_variance_per_site_mpo(*this);
    measurements.energy_variance_per_site_ham   = tools::infinite::measure::energy_variance_per_site_ham(*this);
    measurements.energy_variance_per_site_mom   = tools::infinite::measure::energy_variance_per_site_mom(*this);
    measurements.current_entanglement_entropy   = tools::infinite::measure::current_entanglement_entropy(*this);
}



void  class_infinite_state::swap_AB(){
    MPS->swap_AB();
}
