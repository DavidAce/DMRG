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

std::vector<std::string_view> text::split(std::string_view str, std::string_view dlm) {
    // Here the resulting vector "output" will contain increasingly longer string_views, that are subsets of the path.
    // Note that no string data is allocated here, these are simply views into a string allocated elsewhere.
    // For example if str is "this/is/a/long/path and dlm is "/", the vector will contain the following views
    // [0]: this
    // [1]: is
    // [2]: a
    // [3]: long
    // [4]: path

    // It is very important to note that the resulting views are not null terminated. Therefore, these vector elements
    // **must not** be used as c-style arrays using their .data() member functions.
    std::vector<std::string_view> output;
    size_t                        currentPosition = 0;
    while(currentPosition < str.size()) {
        auto foundpos = str.substr(currentPosition).find_first_of(dlm);
        auto foundstr = str.substr(currentPosition, foundpos);
        output.emplace_back(foundstr);
        if(foundpos == std::string_view::npos) break;
        currentPosition += (foundpos + dlm.size());
    }
    return output;
}

bool text::natcomp(std::string_view sa, std::string_view sb) { return extract_value<long>(sa) < extract_value<long>(sb); }
