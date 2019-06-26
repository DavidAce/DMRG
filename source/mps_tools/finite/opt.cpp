//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <sstream>
#include <spdlog/spdlog.h>
#include <mps_tools/finite/opt.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <general/class_eigsolver.h>
#include <general/arpack_extra/matrix_product_hamiltonian.h>
#include <general/class_svd_wrapper.h>
#include <sim_parameters/nmspc_sim_settings.h>

std::tuple<Eigen::Tensor<std::complex<double>,3>, double>
mpstools::finite::opt::find_optimal_excited_state(const class_finite_chain_state & state, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace,OptType optType){
    mpstools::log->trace("Finding optimal excited state");
    using namespace opt::internals;
    using namespace Textra;
    std::stringstream problem_report;
    auto dims = state.active_dimensions();
    auto size = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());
    problem_report
            << "Starting optimization"
            << std::setprecision(10)
            << "\t mode [ "     << optMode << " ]"
            << "\t space [ "    << optSpace << " ]"
            << "\t type [ "     << optType << " ]"
            << "\t position [ " << state.get_position() << " ]"
            << "\t shape "      << dims << " = [ " << size << " ]" << std::flush;
    mpstools::log->debug(problem_report.str());



    switch (optSpace){
        case OptSpace::FULL:        return internals::subspace_optimization(state, sim_state , optMode, optSpace, optType);
        case OptSpace::PARTIAL:     return internals::subspace_optimization(state, sim_state , optMode, optSpace, optType);
        case OptSpace::DIRECT:      return internals::direct_optimization  (state, sim_state, optType);
    }
}

std::tuple<Eigen::Tensor<std::complex<double>,4>, double> mpstools::finite::opt::find_optimal_ground_state(const class_finite_chain_state & state, const class_simulation_state & sim_state, std::string ritz){
    return internals::ground_state_optimization(state,sim_state,ritz);

}


void mpstools::finite::opt::truncate_theta(Eigen::Tensor<std::complex<double>,3> theta, class_finite_chain_state & state, long chi_, double SVDThreshold){
    mpstools::log->trace("Truncating multitheta");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);

    Eigen::Tensor<class_finite_chain_state::Scalar,4> theta4;
    Eigen::Tensor<class_finite_chain_state::Scalar,3> U = theta;
    Eigen::Tensor<class_finite_chain_state::Scalar,3> V;
    Eigen::Tensor<class_finite_chain_state::Scalar,1> S;
    auto reverse_active_sites = state.active_sites;
    std::reverse(reverse_active_sites.begin(),reverse_active_sites.end());

    while (reverse_active_sites.size() > 1){
        size_t site = reverse_active_sites.front();
        std::cout << "site: "<<  site << std::endl;
        long dim0 = U.dimension(0) / state.get_G(site).dimension(0);
        long dim1 = state.get_G(site).dimension(0);

        theta4 = U
                .reshape(Textra::array4{dim0,dim1, U.dimension(1), U.dimension(2)})
                .shuffle(Textra::array4{0,2,1,3});
        std::tie(U, S, V) = SVD.schmidt(theta4, chi_);
        state.truncation_error[site-1] = SVD.get_truncation_error();
        Eigen::Tensor<std::complex<double>,3> V_L = V.contract(Textra::asDiagonalInversed(state.get_L(site+1)), Textra::idx({2},{0}));
        state.get_G(site) = V_L;
        state.get_L(site) = S;
        reverse_active_sites.pop_front();
    }
    size_t site = reverse_active_sites.front();
    std::cout << "site: "<<  site << std::endl;
    Eigen::Tensor<std::complex<double>,3> L_U = Textra::asDiagonalInversed(state.get_L(site)).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
    state.get_G(site) = L_U;
    mpstools::finite::mps::rebuild_environments(state);
}




void mpstools::finite::opt::truncate_theta(Eigen::Tensor<std::complex<double>,4> theta, class_finite_chain_state & state, long chi_, double SVDThreshold) {
    mpstools::log->trace("Truncating theta");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    auto[U, S, V] = SVD.schmidt(theta, chi_);
    state.truncation_error[state.get_position()+1] = SVD.get_truncation_error();
    state.MPS_C  = S;
    Eigen::Tensor<std::complex<double>,3> L_U = Textra::asDiagonalInversed(state.MPS_L.back().get_L()).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
    Eigen::Tensor<std::complex<double>,3> V_L = V.contract(Textra::asDiagonalInversed(state.MPS_R.front().get_L()), Textra::idx({2},{0}));
    state.MPS_L.back().set_G(L_U);
    state.MPS_R.front().set_G(V_L);
}






void mpstools::finite::opt::internals::reset_timers(){
    t_opt-> reset();
    t_eig-> reset();
    t_ham-> reset();
    t_tot-> reset();
    t_vH2v->reset();
    t_vHv ->reset();
    t_vH2 ->reset();
    t_vH  ->reset();
    t_op  ->reset();
}

template<typename Scalar>
mpstools::finite::opt::internals::MultiComponents<Scalar>::MultiComponents(const class_finite_chain_state & state){
    if constexpr (std::is_same<Scalar,double>::value){
        mpo                          = state.get_multimpo().real();
        mpo2                         = mpo.contract (mpo, Textra::idx({3},{2}));
        auto [envL_cplx,envR_cplx]   = state.get_multienv();
        auto [env2L_cplx,env2R_cplx] = state.get_multienv2();
        envL  = envL_cplx.block.real();         envR  = envR_cplx.block.real();
        env2L = env2L_cplx.block.real();        env2R = env2R_cplx.block.real();
    }

    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        mpo                          = state.get_multimpo();
        mpo2                         = mpo.contract (mpo, Textra::idx({3},{2}));
        auto [envL_cplx,envR_cplx]   = state.get_multienv();
        auto [env2L_cplx,env2R_cplx] = state.get_multienv2();
        envL  = envL_cplx.block;         envR  = envR_cplx.block;
        env2L = env2L_cplx.block;        env2R = env2R_cplx.block;
    }
    dsizes        = state.active_dimensions();
}

template struct mpstools::finite::opt::internals::MultiComponents<double>;
template struct mpstools::finite::opt::internals::MultiComponents<std::complex<double>>;


double mpstools::finite::opt::internals::windowed_func_abs(double x,double window){
    if (std::abs(x) >= window){
        return std::abs(x)-window;
    }else{
        return 0;
    }
}
double mpstools::finite::opt::internals::windowed_grad_abs(double x,double window){
    if (std::abs(x) >= window){
        return sgn(x);
    }else{
        return 0.0;
    }
}



double mpstools::finite::opt::internals::windowed_func_pow(double x,double window){
    if (std::abs(x) >= window){
        return x*x - window*window;
    }else{
        return 0.0;
    }
}
double mpstools::finite::opt::internals::windowed_grad_pow(double x,double window){
    if (std::abs(x) >= window){
        return 2.0*x;
    }else{
        return 0.0;
    }
}


