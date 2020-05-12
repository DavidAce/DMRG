#include "views.h"
#include <model/class_model_finite.h>
#include <state/class_state_finite.h>
#include <state/class_mps_site.h>
#include <model/class_mpo_base.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

using Scalar = std::complex<double>;

const Eigen::Tensor<Scalar, 3> &tools::finite::views::get_multisite_mps(const class_state_finite &state) {
    return tools::finite::views::get_multisite_mps(state, state.active_sites);
}

const Eigen::Tensor<Scalar, 3> &tools::finite::views::get_multisite_mps(const class_state_finite &state, const std::list<size_t> &active_sites) {
    if(state.cache.multisite_mps) return state.cache.multisite_mps.value();
    if(active_sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mps");
    tools::log->trace("Contracting multi theta");
    tools::common::profile::t_mps->tic();
    Eigen::Tensor<Scalar, 3> multisite_mps;
    Eigen::Tensor<Scalar, 3> temp;
    bool                     first = true;
    for(auto &site : active_sites) {
        if(first) {
            multisite_mps = state.get_MPS(site).get_M();
            first         = false;
            continue;
        }
        const auto &M    = state.get_MPS(site).get_M();
        long        dim0 = multisite_mps.dimension(0) * M.dimension(0);
        long        dim1 = multisite_mps.dimension(1);
        long        dim2 = M.dimension(2);
        temp             = multisite_mps.contract(M, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dim0, dim1, dim2});
        multisite_mps    = temp;
    }
    state.cache.multisite_mps = multisite_mps;
    //    state.cache.multisite_mps = temp;
    tools::common::profile::t_mps->tic();
    return state.cache.multisite_mps.value();
}

const Eigen::Tensor<class_model_finite::Scalar, 4> &tools::finite::views::get_multisite_mpo(const class_model_finite &model) {
    return tools::finite::views::get_multisite_mpo(model, model.active_sites);
}

const Eigen::Tensor<class_model_finite::Scalar, 4> &tools::finite::views::get_multisite_mpo(const class_model_finite &model,
                                                                                            const std::list<size_t> & active_sites) {
    if(model.cache.multisite_mpo) return model.cache.multisite_mpo.value();
    if(active_sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo");
    tools::log->trace("Contracting multi mpo");
    tools::common::profile::t_mpo->tic();
    Eigen::Tensor<Scalar, 4> multisite_mpo;
    Eigen::Tensor<Scalar, 4> temp;
    bool                     first = true;
    for(auto &site : active_sites) {
        if(first) {
            multisite_mpo = model.get_MPO(site).MPO();
            first         = false;
            continue;
        }
        const auto &mpo  = model.get_MPO(site).MPO();
        long        dim0 = multisite_mpo.dimension(0);
        long        dim1 = mpo.dimension(1);
        long        dim2 = multisite_mpo.dimension(2) * mpo.dimension(2);
        long        dim3 = multisite_mpo.dimension(3) * mpo.dimension(3);
        temp = multisite_mpo.contract(mpo, Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
        multisite_mpo = temp;
    }
    model.cache.multisite_mpo = multisite_mpo;
    tools::common::profile::t_mpo->toc();
    return model.cache.multisite_mpo.value();
}
