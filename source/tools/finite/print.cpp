//
// Created by david on 2019-01-29.
//

#include <string>
#include <tensors/class_tensors_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/print.h>
using Scalar = std::complex<double>;

void tools::finite::print::dimensions(const class_tensors_finite &tensors) {
    for(size_t pos = 0; pos < tensors.get_length(); pos++) {
        std::string tag;
        if(pos == tensors.get_position()) tag = "<---- Position A";
        if(pos == tensors.get_position() + 1) tag = "<---- Position B";
        const auto &mps  = tensors.state->get_mps_site(pos).dimensions();
        const auto &envl = tensors.edges->get_ene(pos).L.get_block().dimensions();
        const auto &envr = tensors.edges->get_ene(pos).R.get_block().dimensions();
        const auto &mpo  = tensors.model->get_mpo(pos).MPO().dimensions();
        tools::log->info(
            "Pos {:2}: ENVL [{:>3} {:>3} {:>2}] MPS [{:>2} {:>3} {:>3}] MPO [{:>1} {:>1} {:>1} {:>1}] ENVR [{:>3} {:>3} {:>2}] {}", pos,
            envl[0], envl[1], envl[2], mps[0], mps[1], mps[2], mpo[0], mpo[1], mpo[2], mpo[3], envr[0],
            envr[1], envr[2], tag);
        if(tensors.state->get_mps_site(pos).isCenter())
            tools::log->info("Pos {:2}: L [{:^4}] {:>70}", pos, tensors.state->get_mps_site(pos).get_L().dimension(0), "<---- Center");
    }
}

void tools::finite::print::model(const class_model_finite &model) {
    model.get_mpo(0).print_parameter_names();
    for(size_t pos = 0; pos < model.get_length(); pos++) model.get_mpo(pos).print_parameter_values();
}
