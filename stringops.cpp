#include "stringops.hpp"
namespace stringops {

bool endswith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<string> split(const string& s, const char* delims, int nsplit) {
    std::vector<string> outp;

    if (s.empty()) { return outp; }

    std::vector<int> delimsites;
    for (size_t i = 0; i < s.size(); ++i) {
        if (nsplit > -1 && delimsites.size() == nsplit) break;
        for (const char* d = delims; *d != '\0'; ++d) {
            if (s[i] == *d) delimsites.push_back(i);

        }
    }


    size_t last_delim = 0;
    for (int i : delimsites) {
        string tok = s.substr(last_delim, i - last_delim);
        outp.push_back(tok);
        last_delim = i + 1;
    }
    string tok = s.substr(last_delim);
    outp.push_back(tok);
    return outp;
}

string join(const std::vector<string>& toks, const string& delim) {
    size_t total_size = delim.size() * (toks.size() - 1);
    for (auto s : toks) { total_size += s.size(); }

    string outp;
    outp.reserve(total_size);

    
    for (size_t i=0; i < toks.size(); i++) {
        string s = toks[i];
        outp.append(s);
        if (i != (toks.size()-1)) outp.append(delim);
    }
    return outp;
}

string get_single_token(const string& s, const char sep, size_t tokidx) {
    size_t s_size = s.length();

    size_t startidx = 0;
    size_t curtokidx = 0;
    while (curtokidx < tokidx) {
        if (startidx == s.size()) { throw std::out_of_range("not enough tokens"); }
        if (s[startidx] == sep) {
            curtokidx++;
            if (curtokidx == tokidx) {startidx++; break;}
        }
        startidx++;
    }

    size_t stopidx = startidx + 1;
    while (s[stopidx] != sep && stopidx != s_size) {
        stopidx++;
    }
    return s.substr(startidx, stopidx - startidx);
}

bool startswith(const string& s, const string& start) {
    if (s.empty() || start.empty()) return false;
    else if (start.size() > s.size()) return false;
    else if (start.size() == s.size()) return (s.compare(start) == 0);

    for (size_t i = 0; i < start.size(); ++i) {
        if (s[i] != start[i]) return false;
    }
    return true;
}

}