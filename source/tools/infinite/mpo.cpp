//
// Created by david on 2019-10-24.
//

#include <tensors/model/class_mpo_base.h>
#include <tensors/model/class_mpo_factory.h>
#include <tensors/model/class_model_infinite.h>
#include <tools/common/log.h>
#include <tools/infinite/mpo.h>

void tools::infinite::mpo::initialize(class_model_infinite & model, ModelType model_type){
    tools::log->trace("Initializing mpo");
    //Generate MPO
    model.model_type = model_type;
    model.HA = class_mpo_factory::create_mpo(0,model_type);
    model.HB = class_mpo_factory::create_mpo(1,model_type);
}

void tools::infinite::mpo::randomize(class_model_infinite &model) {
    tools::log->trace("Setting random fields in MPO's");
    model.HA->randomize_hamiltonian();
    model.HB->randomize_hamiltonian();
    std::vector<class_mpo_base::TableMap> all_params;

    all_params.push_back(model.HA->get_parameters());
    all_params.push_back(model.HB->get_parameters());

    model.HA->set_averages(all_params);
    model.HB->set_averages(all_params);

}