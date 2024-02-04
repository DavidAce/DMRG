#include "tid.h"
namespace tid {
    token::token(ur &t_) : t(t_) {
        if(t.lvl == level::disabled) return;
        t.tic();
    }
    token::token(ur &t_, std::string_view prefix_) : t(t_) {
        if(t.lvl == level::disabled) return;
        temp_prefix = prefix_;
        tid::internal::ur_prefix_push_back(temp_prefix);
        t.tic();
    }

    token::~token() noexcept {
        if(t.lvl == level::disabled) return;
        try {
            if(t.is_measuring) t.toc();
            tid::internal::ur_prefix_pop_back(temp_prefix);
        } catch(const std::exception &ex) { fprintf(stderr, "tid: error in token destructor for tid::ur [%s]: %s", t.get_label().c_str(), ex.what()); }
    }

    void token::tic() noexcept {
        if(t.lvl == level::disabled) return;
        t.tic();
    }
    void token::toc() noexcept {
        if(t.lvl == level::disabled) return;
        t.toc();
        tid::internal::ur_prefix_pop_back(temp_prefix);
    }
    ur &token::ref() noexcept { return t; }
    ur *token::operator->() noexcept { return &t; }

}
