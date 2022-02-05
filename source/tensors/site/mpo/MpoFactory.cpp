#include "MpoFactory.h"
#include "IsingSelfDual.h"
#include "IsingMajorana.h"
#include "IsingRandomField.h"
#include "LBit.h"
#include "MpoSite.h"

std::unique_ptr<MpoSite> MpoFactory::create_mpo(size_t position, ModelType model_type) {
    switch(model_type) {
        case ModelType::ising_tf_rf: return std::make_unique<IsingRandomField>(model_type, position);
        case ModelType::ising_sdual: return std::make_unique<IsingSelfDual>(model_type, position);
        case ModelType::ising_majorana: return std::make_unique<IsingMajorana>(model_type, position);
        case ModelType::lbit: return std::make_unique<LBit>(model_type, position);
        default: throw std::runtime_error(fmt::format("Wrong model type: [{}]", enum2sv(model_type)));
    }
}

std::unique_ptr<MpoSite> MpoFactory::clone(std::unique_ptr<MpoSite> other) { return other->clone(); }
