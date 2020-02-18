//
// Created by david on 2018-07-04.
//

#include "class_model_factory.h"
#include "class_model_base.h"
#include "class_model_parameters.h"
#include "class_selfdual_tf_rf_ising.h"
#include "class_tf_ising.h"
#include <math/nmspc_math.h>

std::unique_ptr<class_model_base> class_model_factory::create_mpo(size_t position, std::string_view model_type_str) {
    if(model_type_str == std::string("tf_ising")) return std::make_unique<class_tf_ising>(position);
    if(model_type_str == std::string("selfdual_tf_rf_ising"))
        return std::make_unique<class_selfdual_tf_rf_ising>(position);
    else
        throw std::runtime_error("Wrong model: [ " + std::string(model_type_str) + " ]");
}

std::unique_ptr<class_model_base> class_model_factory::clone(std::unique_ptr<class_model_base> other) { return other->clone(); }
