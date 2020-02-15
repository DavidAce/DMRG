//
// Created by david on 2019-06-24.
//

#include <tools/infinite/debug.h>
#include <tools/infinite/print.h>
#include <tools/infinite/measure.h>
#include <tools/common/log.h>
#include <state/class_state_infinite.h>

void tools::infinite::debug::check_integrity(const class_state_infinite & state)
{
    tools::log->debug("Checking integrity...");
   try{
        state.assert_positions();
        check_integrity_of_mps(state);
    }catch(std::exception & ex){
        tools::infinite::print::print_state(state) ;
        throw std::runtime_error("Integrity check of state failed: " + std::string(ex.what()));
    }
}




void tools::infinite::debug::check_integrity_of_mps(const class_state_infinite &state){
    tools::log->debug("Checking integrity of MPS...");
    tools::log->debug("Checking norms");
    auto norm_block = tools::infinite::measure::norm(state);
    if(std::abs(norm_block - 1.0) > 1e-10) {
        tools::log->warn("Norm of state far from unity: {}", norm_block);
        throw std::runtime_error("Norm of state too far from unity: " + std::to_string(norm_block));
    }
    tools::log->debug("MPS OK");
}



