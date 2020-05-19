//
// Created by david on 2018-07-04.
//

#pragma once

#include "class_mpo_site.h"
#include <iostream>
#include <memory>

class class_mpo_factory {
    public:
    static std::unique_ptr<class_mpo_site> create_mpo(size_t position, ModelType model_type);

    template<typename... T>
    static std::unique_ptr<class_mpo_site> create_mpo(size_t position, ModelType model_type, T... args) {
        auto mpo = create_mpo(position, model_type);
        mpo->set_parameters(args...);
        return mpo;
    }
    //    static std::unique_ptr<class_hamiltonian_h5log_base> create_table(std::string model_type_str);
    static std::unique_ptr<class_mpo_site> clone(std::unique_ptr<class_mpo_site> other);
};
