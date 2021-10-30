#include "prof.h"
#include <config/settings.h>
#include <tid/tid.h>
#include <tools/common/log.h>
void tools::common::timer::print_timers() {
    if(settings::timer::on) {
        tid::level lvl = settings::timer::level;

        if(settings::itebd::on)
            for(const auto &t : tid::get_tree("iTEBD"))
                if(t->get_level() <= lvl) tools::log->info("{}", t.str());

        if(settings::flbit::on)
            for(const auto &t : tid::get_tree("fLBIT"))
                if(t->get_level() <= lvl) tools::log->info("{}", t.str());

        if(settings::idmrg::on)
            for(const auto &t : tid::get_tree("iDMRG"))
                if(t->get_level() <= lvl) tools::log->info("{}", t.str());

        if(settings::fdmrg::on)
            for(const auto &t : tid::get_tree("fDMRG"))
                if(t->get_level() <= lvl) tools::log->info("{}", t.str());

        if(settings::xdmrg::on) {
            for(const auto &t : tid::get_tree("fDMRG"))
                if(t->get_level() <= lvl) tools::log->info("{}", t.str());
            for(const auto &t : tid::get_tree("xDMRG"))
                if(t->get_level() <= lvl) tools::log->info("{}", t.str());
        }
    }
}