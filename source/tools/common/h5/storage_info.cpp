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

StoragePolicy StorageInfo::get_state_storage_policy() const{
    if(state_name == "state_init") return settings::storage::mps::state_init::policy;
    if(state_name == "state_emid") return settings::storage::mps::state_emid::policy;
    if(state_name == "state_emin") return settings::storage::mps::state_emin::policy;
    if(state_name == "state_emax") return settings::storage::mps::state_emax::policy;
    if(state_name == "state_real") return settings::storage::mps::state_real::policy;
    if(state_name == "state_lbit") return settings::storage::mps::state_lbit::policy;

    // xDMRG states are numbered like "state_#". We should therefore get the policy simply based on the algorithm
    if(algo_type == AlgorithmType::xDMRG) return settings::storage::mps::state_emid::policy;
    throw except::logic_error("Failed to get the storage policy for state name: [{}]", state_name);
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
        case StorageEvent::ITERATION: return fmt::format("{}/mps/iter_{}", state_prefix, iter);
        case StorageEvent::INIT: return fmt::format("{}/mps", state_prefix);
        case StorageEvent::FINISHED: return fmt::format("{}/mps", state_prefix);
        case StorageEvent::PROJECTION: return fmt::format("{}/mps/proj", state_prefix);
        case StorageEvent::BOND_UPDATE: return fmt::format("{}/mps/bond_{}", state_prefix, bond_lim);
        case StorageEvent::TRNC_UPDATE: return fmt::format("{}/mps/trnc_{:.2e}", state_prefix, trnc_lim);
        case StorageEvent::FES_STEP: return fmt::format("{}/mps/fes_{}", state_prefix, bond_lim);
        case StorageEvent::EMIN_STATE: return fmt::format("{}/mps", state_prefix);
        case StorageEvent::EMAX_STATE: return fmt::format("{}/mps", state_prefix);
        case StorageEvent::MODEL: throw except::logic_error("get_mps_prefix(): Invalid event: [MODEL]");
        case StorageEvent::NONE: throw except::logic_error("get_mps_prefix(): Invalid event: [NONE]");
        default: throw except::logic_error("Unhandled event: [{}]", enum2sv(storage_event));
    }
}

StorageInfo::StorageInfo(const AlgorithmStatus &status, std::string_view state_name, StorageEvent event)
    : iter(status.iter),                                     //
      step(status.step), position(status.position),          //
      direction(status.direction),                           //
      bond_lim(status.bond_lim),                             //
      bond_max(status.bond_max),                             //
      trnc_lim(status.trnc_lim),                             //
      algo_type(status.algo_type),                           //
      algo_name(status.algo_type_sv()),                      //
      storage_event(status.event),                           //
      state_name(state_name),                                //
      algorithm_has_finished(status.algorithm_has_finished), //
      algorithm_has_succeeded(status.algorithm_has_succeeded) {
    if(storage_event == StorageEvent::NONE) storage_event = event;
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
