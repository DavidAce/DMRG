//
// Created by david on 2019-02-20.
//

#include "print.h"
#include <tensors/state/class_state_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/edges/class_env_ene.h>

void tools::infinite::print::print_hamiltonians(const class_model_infinite &model) {
    model.get_mpo_siteA().print_parameter_names();
    model.get_mpo_siteA().print_parameter_values();
    model.get_mpo_siteB().print_parameter_values();
}

void tools::infinite::print::dimensions(const class_tensors_infinite &tensors) {
    const auto &mpsa = tensors.state->get_mps_siteA().dimensions();
    const auto &mpsb = tensors.state->get_mps_siteB().dimensions();
    std::string strA = "Site A ";
    std::string strB = "Site B ";
    if(tensors.state){
        strA += fmt::format("MPS [{:>2} {:>3} {:>3}] ",mpsa[0], mpsa[1], mpsa[2]);
        strB += fmt::format("MPS [{:>2} {:>3} {:>3}] ",mpsb[0], mpsb[1], mpsb[2]);
    }

    if(tensors.edges){
        const auto &envl = tensors.edges->get_ene().L.get_block().dimensions();
        const auto &envr = tensors.edges->get_ene().R.get_block().dimensions();
        strA = fmt::format("ENV [{:>3} {:>3} {:>2}] ",envl[0], envl[1], envl[2]);
        strB = fmt::format("ENV [{:>3} {:>3} {:>2}] ",envr[0], envr[1], envr[2]);
    }
    if(tensors.model){
        const auto &mpoa  = tensors.model->get_mpo_siteA().MPO().dimensions();
        const auto &mpob  = tensors.model->get_mpo_siteB().MPO().dimensions();
        strA = fmt::format("MPO [{:>1} {:>1} {:>1} {:>1}] ",mpoa[0], mpoa[1], mpoa[2], mpoa[3]);
        strB = fmt::format("MPO [{:>1} {:>1} {:>1} {:>1}] ",mpob[0], mpob[1], mpob[2], mpob[3]);
    }
    tools::log->info(strA);
    tools::log->info(strB);
}


