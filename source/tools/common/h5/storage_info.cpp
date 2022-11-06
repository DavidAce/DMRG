#include "storage_info.h"
#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "fmt/core.h"
#include <string>
void StorageInfo::assert_well_defined() const {
    std::string msg;
    if(iter == -1ul) msg.append(" iter\n");
    if(step == -1ul) msg.append(" step\n");
    if(position == -1) msg.append(" position\n");
    if(direction == 0) msg.append(" direction\n");
    if(bond_lim == -1) msg.append(" bond_lim\n");
    if(bond_max == -1) msg.append(" bond_max\n");
    if(trnc_lim == -1) msg.append(" trnc_lim\n");
    if(algo_type == AlgorithmType::ANY) msg.append(" algo_type\n");
    if(algo_name.empty()) msg.append(" algo_name\n");
    if(state_name.empty()) msg.append(" state_name\n");
    if(not msg.empty()) throw except::logic_error("StorageInfo is not well defined. Missing fields: \n{}", msg);
}

std::string StorageInfo::get_state_prefix() const {
    assert_well_defined();
    return fmt::format("{}/{}", algo_name, state_name);
}
std::string StorageInfo::get_mps_prefix() const {
    assert_well_defined();
    std::string tag;
    auto        state_prefix = get_state_prefix();
    switch(storage_event) {
        case StorageEvent::ITER_STATE: return fmt::format("{}/mps/iter_{}", state_prefix, iter);
        case StorageEvent::INIT_STATE: return fmt::format("{}/mps/init", state_prefix);
        case StorageEvent::LAST_STATE: return fmt::format("{}/mps", state_prefix);
        case StorageEvent::PROJ_STATE: return fmt::format("{}/mps/proj", state_prefix);
        case StorageEvent::BOND_INCREASE: return fmt::format("{}/{}/mps/bond_{}", state_prefix, bond_lim);
        case StorageEvent::TRNC_DECREASE: return fmt::format("{}/{}/mps/trnc_{:.2e}", state_prefix, trnc_lim);
        case StorageEvent::FES_STATE: return fmt::format("{}/{}/mps/fes_{}", state_prefix, bond_lim);
        case StorageEvent::EMIN_STATE: return fmt::format("{}/state_emin/mps", algo_name);
        case StorageEvent::EMAX_STATE: return fmt::format("{}/state_emax/mps", algo_name);
        case StorageEvent::MODEL: throw except::logic_error("get_mps_prefix(): Invalid event: [MODEL]");
        case StorageEvent::NONE: throw except::logic_error("get_mps_prefix(): Invalid event: [NONE]");
        default: throw except::logic_error("Unhandled event: [{}]", enum2sv(storage_event));
    }
}

StorageInfo::StorageInfo(const AlgorithmStatus &status, std::string_view state_name, StorageEvent event)
    : iter(status.iter), step(status.step), position(status.position), direction(status.direction), bond_lim(status.bond_lim), bond_max(status.bond_max),
      trnc_lim(status.trnc_lim), algo_type(status.algo_type), algo_name(status.algo_type_sv()), storage_event(status.event), state_name(state_name) {
    if(storage_event == StorageEvent::NONE) storage_event = event;

    switch(storage_event) {
        case StorageEvent::ITER_STATE: {
            if(settings::storage::storage_interval == 0) break;
            if(status.iter % settings::storage::storage_interval != 0) break;
            storage_level = settings::storage::storage_level_iter_state;
            break;
        }
        case StorageEvent::INIT_STATE: storage_level = settings::storage::storage_level_init_state; break;
        case StorageEvent::LAST_STATE: storage_level = settings::storage::storage_level_last_state; break;
        case StorageEvent::EMIN_STATE: storage_level = settings::storage::storage_level_emin_state; break;
        case StorageEvent::EMAX_STATE: storage_level = settings::storage::storage_level_emax_state; break;
        case StorageEvent::PROJ_STATE: storage_level = settings::storage::storage_level_proj_state; break;
        case StorageEvent::BOND_INCREASE: storage_level = settings::storage::storage_level_bond_state; break;
        case StorageEvent::TRNC_DECREASE: storage_level = settings::storage::storage_level_trnc_state; break;
        case StorageEvent::FES_STATE: storage_level = settings::storage::storage_level_fes_state; break;
        case StorageEvent::MODEL: storage_level = settings::storage::storage_level_model; break;
        case StorageEvent::NONE: break;
    }
    assert_well_defined();
}

StorageInfo::~StorageInfo() noexcept { storage_event = StorageEvent::NONE; }

StorageAttrs::StorageAttrs(const StorageInfo &sinfo) {
    iter          = sinfo.iter;
    step          = sinfo.step;
    bond_lim      = sinfo.bond_lim;
    bond_max      = sinfo.bond_max;
    trnc_lim      = sinfo.trnc_lim;
    storage_event = sinfo.storage_event;
    storage_level = sinfo.storage_level;
}
bool StorageAttrs::operator==(const StorageAttrs &sinfo) {
    return iter == sinfo.iter and step == sinfo.step and bond_lim == sinfo.bond_lim and bond_max == sinfo.bond_max and trnc_lim == sinfo.trnc_lim and
           storage_event == sinfo.storage_event and storage_level == sinfo.storage_level;
}

bool StorageAttrs::operator==(const StorageInfo &sinfo) {
    return iter == sinfo.iter and step == sinfo.step and bond_lim == sinfo.bond_lim and bond_max == sinfo.bond_max and trnc_lim == sinfo.trnc_lim and
           storage_event == sinfo.storage_event and storage_level == sinfo.storage_level;
}
