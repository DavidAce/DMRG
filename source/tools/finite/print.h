#pragma once

class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
namespace tools::finite::print {
    extern void dimensions(const StateFinite &state);
    extern void dimensions(const EdgesFinite &edges);
    extern void dimensions(const TensorsFinite &tensors);
    extern void model(const ModelFinite &model);

}