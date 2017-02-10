#ifndef STRINGOPS_HPP
#define STRINGOPS_HPP

#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <memory>



namespace stringops {
    using std::string;
    bool endswith(const string& s, const string& suffix);
    string join(const std::vector<string>& v, const string& delim);
    std::vector<string> split(const string& s, const char* delim, int nsplit=-1);
    bool startswith(const string& s, const string& start);
    string get_single_token(const string& s, const char sep, size_t tokidx);
}


#endif