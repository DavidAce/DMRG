#include "tid.h"
namespace tid {
    token::token(ur &t_) : t(t_) { t.tic(); }
    token::token(ur &t_, std::string_view prefix_) : t(t_) {
        if(tid::internal::ur_prefix.empty()) {
            temp_prefix = std::string(prefix_);
        } else {
            temp_prefix = "." + std::string(prefix_);
        }
        tid::internal::ur_prefix.append(temp_prefix);
        t.tic();
    }

    token::~token() noexcept {
        try {
            if(t.is_measuring) t.toc();
            if(not temp_prefix.empty()) {
                auto pos = tid::internal::ur_prefix.rfind(temp_prefix);
                if(pos != std::string::npos) { tid::internal::ur_prefix.erase(pos, temp_prefix.size()); }
            }
        } catch(const std::exception &ex) { fprintf(stderr, "tid: error in token destructor for tid::ur [%s]: %s", t.get_label().c_str(), ex.what()); }
    }

    void token::tic() noexcept { t.tic(); }
    void token::toc() noexcept {
        t.toc();
        if(not temp_prefix.empty()) {
            auto pos = tid::internal::ur_prefix.rfind(temp_prefix);
            if(pos != std::string::npos) { tid::internal::ur_prefix.erase(pos, temp_prefix.size()); }
        }
    }
    ur &token::ref() noexcept { return t; }
    ur *token::operator->() noexcept { return &t; }

}
