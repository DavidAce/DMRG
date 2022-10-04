#include "text.h"
bool text::endsWith(std::string_view str, std::string_view suffix) {
    return str.size() >= suffix.size() and 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

bool text::startsWith(std::string_view str, std::string_view prefix) { return str.size() >= prefix.size() and 0 == str.compare(0, prefix.size(), prefix); }

std::string text::replace(std::string_view str, std::string_view from, std::string_view to) {
    std::string res(str);
    if(from.empty()) return res;
    size_t start_pos = 0;
    while((start_pos = res.find(from, start_pos)) != std::string::npos) {
        res.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
    return res;
}

bool text::natcomp(std::string_view sa, std::string_view sb) { return extract_value<long>(sa) < extract_value<long>(sb); }
