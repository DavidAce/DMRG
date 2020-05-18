#include "views.h"
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/model/class_mpo_base.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

using Scalar = std::complex<double>;

const Eigen::Tensor<Scalar, 3> &tools::finite::views::get_multisite_tensor(const class_state_finite &state, std::optional<std::list<size_t>> active_sites) {
    if(not active_sites) active_sites = state.active_sites;
    if(active_sites.value().empty()) throw std::runtime_error("No active sites on which to build a multisite mps");
    if(active_sites.value() == state.active_sites and state.cache.multisite_tensor) return state.cache.multisite_tensor.value();
    tools::log->trace("Contracting multi theta");
    tools::common::profile::t_mps->tic();
    Eigen::Tensor<Scalar, 3> multisite_tensor;
    Eigen::Tensor<Scalar, 3> temp;
    bool                     first = true;
    for(auto &site : active_sites.value()) {
        if(first) {
            multisite_tensor = state.get_mps_site(site).get_M();
            first         = false;
            continue;
        }
        const auto &M    = state.get_mps_site(site).get_M();
        long        dim0 = multisite_tensor.dimension(0) * M.dimension(0);
        long        dim1 = multisite_tensor.dimension(1);
        long        dim2 = M.dimension(2);
        temp             = multisite_tensor.contract(M, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dim0, dim1, dim2});
        multisite_tensor    = temp;
    }
    state.cache.multisite_tensor = multisite_tensor;
    //    state.cache.multisite_tensor = temp;
    tools::common::profile::t_mps->tic();
    return state.cache.multisite_tensor.value();
}


const Eigen::Tensor<class_model_finite::Scalar, 4> &tools::finite::views::get_multisite_tensor(const class_model_finite &model,
                                                                                            std::optional<std::list<size_t>> active_sites) {
    if(not active_sites) active_sites = model.active_sites;
    if(active_sites.value().empty()) throw std::runtime_error("No active sites on which to build a multisite mpo");
    if(active_sites.value() == model.active_sites and model.cache.multisite_tensor) return model.cache.multisite_tensor.value();
    tools::log->trace("Contracting multi mpo");
    tools::common::profile::t_mpo->tic();
    Eigen::Tensor<Scalar, 4> multisite_tensor;
    Eigen::Tensor<Scalar, 4> temp;
    bool                     first = true;
    for(auto &site : active_sites.value()) {
        if(first) {
            multisite_tensor = model.get_mpo(site).MPO();
            first         = false;
            continue;
        }
        const auto &site_tensor = model.get_mpo(site).MPO();
        long        dim0 = multisite_tensor.dimension(0);
        long        dim1 = site_tensor.dimension(1);
        long        dim2 = multisite_tensor.dimension(2) * site_tensor.dimension(2);
        long        dim3 = multisite_tensor.dimension(3) * site_tensor.dimension(3);
        temp = multisite_tensor.contract(site_tensor, Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
        multisite_tensor = temp;
    }
    model.cache.multisite_tensor = multisite_tensor;
    tools::common::profile::t_mpo->toc();
    return model.cache.multisite_tensor.value();
}
