#pragma once

#include <memory>
#include <string_view>

class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_algorithm_status;
namespace h5pp {
    class File;
}

namespace tools::finite{
    class mwrap {
        // Mutable wrapper
        public:
        class_state_finite &     state;
        class_model_finite &     model;
        class_edges_finite &     edges;
        class_algorithm_status &status;
        h5pp::File &             file;
        mwrap(class_state_finite &state,  class_model_finite &model, class_edges_finite &edges, class_algorithm_status &status,
                    h5pp::File &file) : state(state), model(model), edges(edges), status(status), file(file){};
    };

    class cwrap {
        // Const wrapper
        public:
        const class_state_finite &     state;
        const class_model_finite &     model;
        const class_edges_finite &     edges;
        const class_algorithm_status &status;
        const h5pp::File &             file;
        cwrap(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges, const class_algorithm_status &status,
                     const h5pp::File &file)
            : state(state), model(model), edges(edges), status(status), file(file){};
    };


    class wwrap {
        // Writeable wrapper
        public:
        const class_state_finite &     state;
        const class_model_finite &     model;
        const class_edges_finite &     edges;
        const class_algorithm_status &status;
        h5pp::File &                   file;
        wwrap(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges, const class_algorithm_status &status,
                     h5pp::File &file)
            : state(state), model(model), edges(edges), status(status), file(file){};
    };


}

