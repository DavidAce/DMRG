//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <sstream>
#include <spdlog/spdlog.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <general/class_tic_toc.h>
#include <LBFGS.h>




std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::find_optimal_excited_state(const class_superblock & superblock, class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace){
    MPS_Tools::log->trace("Finding optimal excited state");
    using namespace Opt::internals;
    std::stringstream problem_report;
    auto dims = superblock.dimensions();
    problem_report
            << "Starting optimization"
            << std::setprecision(10)
            << "\t mode [ "    << optSpace << " ]"
            << "\t position [ "    << superblock.get_position() << " ]"
            << "\t chi [ "    << superblock.get_chi() << " ]"
            << "\t shape [ "    << dims[0] << " " << dims[1] << " " << dims[2]<< " " <<dims[3]  << " ] = [ " << dims[0]*dims[1]*dims[2]*dims[3] << " ]" << std::flush;
    MPS_Tools::log->debug(problem_report.str());



    switch (optSpace){
        case OptSpace::DIRECT:      return internals::direct_optimization(superblock,sim_state);
        case OptSpace::GUIDED:      return internals::guided_optimization(superblock,sim_state);
        case OptSpace::CPPOPTLIB:   return internals::cppoptlib_optimization(superblock,sim_state);
        case OptSpace::FULL:        return internals::subspace_optimization(superblock, sim_state , optMode, optSpace);
        case OptSpace::PARTIAL:     return internals::subspace_optimization(superblock, sim_state , optMode, optSpace);
    }
}

namespace MPS_Tools::Finite::Opt::internals{
    std::shared_ptr<LBFGSpp::LBFGSParam<double>> params = std::make_shared<LBFGSpp::LBFGSParam<double>>();
}



void MPS_Tools::Finite::Opt::internals::initialize_params(){
    using namespace LBFGSpp;
    // READ HERE http://pages.mtu.edu/~msgocken/ma5630spring2003/lectures/lines/lines/node3.html
    // I think c1 corresponds to ftol, and c2 corresponds to wolfe
    params->max_iterations = 1000;
    params->max_linesearch = 60; // Default is 20. 5 is really bad, 80 seems better.
    params->m              = 8;     // Default is 6
    params->past           = 1;     // Or perhaps it was this one that helped.
    params->epsilon        = 1e-2;  // Default is 1e-5.
    params->delta          = 1e-6; // Default is 0. Trying this one instead of ftol.
    params->ftol           = 1e-4;  // Default is 1e-4. this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
    params->wolfe          = 0.90;   // Default is 0.9
    params->min_step       = 1e-40;
    params->max_step       = 1e+40;
    params->linesearch     = LINE_SEARCH_ALGORITHM::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
}




namespace MPS_Tools::Finite::Opt::internals{
    std::shared_ptr<class_tic_toc>  t_opt  =  std::make_shared<class_tic_toc>(true,5,"t_opt ");
    std::shared_ptr<class_tic_toc>  t_eig  =  std::make_shared<class_tic_toc>(true,5,"t_eig ");
    std::shared_ptr<class_tic_toc>  t_ham  =  std::make_shared<class_tic_toc>(true,5,"t_ham ");
    std::shared_ptr<class_tic_toc>  t_tot  =  std::make_shared<class_tic_toc>(true,5,"t_tot ");
    std::shared_ptr<class_tic_toc>  t_vH2v =  std::make_shared<class_tic_toc>(true,5,"t_vH2v");
    std::shared_ptr<class_tic_toc>  t_vHv  =  std::make_shared<class_tic_toc>(true,5,"t_vHv ");
    std::shared_ptr<class_tic_toc>  t_vH2  =  std::make_shared<class_tic_toc>(true,5,"t_vH2 ");
    std::shared_ptr<class_tic_toc>  t_vH   =  std::make_shared<class_tic_toc>(true,5,"t_vH  ");
    std::shared_ptr<class_tic_toc>  t_op   =  std::make_shared<class_tic_toc>(true,5,"t_op  ");
}



void MPS_Tools::Finite::Opt::internals::reset_timers(){
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




std::pair<Eigen::VectorXd,double> MPS_Tools::Finite::Opt::internals::get_vH_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents){
    auto vH = MPS_Tools::Finite::Opt::internals::get_vH(v,superComponents);
    t_vHv->tic();
    auto vHv = vH.dot(v);
    t_vHv->toc();
    return std::make_pair(vH,vHv);
}

std::pair<Eigen::VectorXd,double> MPS_Tools::Finite::Opt::internals::get_vH2_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents){
    auto vH2 = MPS_Tools::Finite::Opt::internals::get_vH2(v,superComponents);
    t_vH2v->tic();
    auto vH2v = vH2.dot(v);
    t_vH2v->toc();
    return std::make_pair(vH2,vH2v);
}


Eigen::VectorXd MPS_Tools::Finite::Opt::internals::get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superComponents.dsizes);
    t_vH2->tic();
    Eigen::Tensor<double, 4> vH2 =
            superComponents.Lblock2
                    .contract(theta,                            Textra::idx({0},{1}))
                    .contract(superComponents.HAHB2,            Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(superComponents.Rblock2,          Textra::idx({1,2,4},{0,2,3}))
                    .shuffle(Textra::array4{1,0,2,3});
    t_vH2->toc();
////
//    t_vH2->tic();
//    Eigen::Tensor<double, 4> vH2 =
//            theta
//            .contract(superComponents.Lblock2HAHA,Textra::idx({1,0},{0,1}))
//            .contract(superComponents.Rblock2HBHB,Textra::idx({1,0,2,3},{0,1,2,3}));
//    t_vH2->toc();


    return Eigen::Map<Eigen::VectorXd>(vH2.data(),vH2.size());
}

Eigen::VectorXd MPS_Tools::Finite::Opt::internals::get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superComponents.dsizes);
    t_vH->tic();
    Eigen::Tensor<double, 4> vH =
            superComponents.Lblock
                    .contract(theta,                               Textra::idx({0},{1}))
                    .contract(superComponents.HAHB,                Textra::idx({1,2,3},{0,1,4}))
                    .contract(superComponents.Rblock ,             Textra::idx({1,3},{0,2}))
                    .shuffle(Textra::array4{1,0,2,3});
    t_vH->toc();

    return Eigen::Map<Eigen::VectorXd>(vH.data(),vH.size());
}
