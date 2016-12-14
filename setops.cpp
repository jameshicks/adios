#include "setops.hpp"


namespace setops {

    vector<int> union_(const vector<int>& a, const vector<int>& b) {
        vector<int> outp; 
        std::set_union(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    vector<int> intersection(const vector<int>& a, const vector<int>& b) {
        vector<int> outp;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    vector<int> symmetric_difference(const vector<int>& a, const vector<int>& b) {
        vector<int> outp;
        std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    vector<int> difference(const vector<int>& a, const vector<int>& b) {
        vector<int> outp;
        std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    vector<int> multi_union(std::initializer_list<vector<int>> inps) {
        vector<int> outp;
        for (auto it = inps.begin(); it != inps.end(); ++it) {
            outp = union_(outp, *it);
        }
        return outp;
    }

    vector<int> multi_intersection(std::initializer_list<vector<int>> inps) {
        vector<int> outp;
        auto it = inps.begin();
        outp = *it;

        for (; it != inps.end(); ++it) {
            outp = intersection(outp, *it);
        }
        return outp;
    }



    set<int> union_(const set<int>& a, const set<int>& b) {
        set<int> outp; 
        std::set_union(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    set<int> intersection(const set<int>& a, const set<int>& b) {
        set<int> outp;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    set<int> symmetric_difference(const set<int>& a, const set<int>& b) {
        set<int> outp;
        std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    set<int> difference(const set<int>& a, const set<int>& b) {
        set<int> outp;
        std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    set<int> multi_union(std::initializer_list<set<int>> inps) {
        set<int> outp;
        for (auto it = inps.begin(); it != inps.end(); ++it) {
            outp = union_(outp, *it);
        }
        return outp;
    }

    set<int> multi_intersection(std::initializer_list<set<int>> inps) {
        set<int> outp;
        auto it = inps.begin();
        outp = *it;

        for (; it != inps.end(); ++it) {
            outp = intersection(outp, *it);
        }
        return outp;
    }

    unordered_set<int> union_(const unordered_set<int>& u_a, const unordered_set<int>& u_b) {
        // unordered_sets can't be used in set operations, so we have to 
        // convert them to ordered sets before performing operations.
        set<int> a;
        for (auto v : u_a) { a.insert(v); }
        set<int> b;
        for (auto v : u_b) { b.insert(v); }

        unordered_set<int> outp; 
        std::set_union(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    unordered_set<int> intersection(const unordered_set<int>& u_a, const unordered_set<int>& u_b) {
        // unordered_sets can't be used in set operations, so we have to 
        // convert them to ordered sets before performing operations.
        set<int> a;
        for (auto v : u_a) { a.insert(v); }
        set<int> b;
        for (auto v : u_b) { b.insert(v); }

        unordered_set<int> outp;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    unordered_set<int> symmetric_difference(const unordered_set<int>& u_a, const unordered_set<int>& u_b) {
        // unordered_sets can't be used in set operations, so we have to 
        // convert them to ordered sets before performing operations.
        set<int> a;
        for (auto v : u_a) { a.insert(v); }
        set<int> b;
        for (auto v : u_b) { b.insert(v); }

        unordered_set<int> outp;
        std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    unordered_set<int> difference(const unordered_set<int>& u_a, const unordered_set<int>& u_b) {
        // unordered_sets can't be used in set operations, so we have to 
        // convert them to ordered sets before performing operations.
        set<int> a;
        for (auto v : u_a) { a.insert(v); }
        set<int> b;
        for (auto v : u_b) { b.insert(v); }

        unordered_set<int> outp;
        std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(outp, outp.begin()));
        return outp;
    }

    unordered_set<int> multi_union(std::initializer_list<unordered_set<int>> inps) {
        unordered_set<int> outp;
        for (auto it = inps.begin(); it != inps.end(); ++it) {
            outp = union_(outp, *it);
        }
        return outp;
    }

    unordered_set<int> multi_intersection(std::initializer_list<unordered_set<int>> inps) {
        unordered_set<int> outp;
        auto it = inps.begin();
        outp = *it;

        for (; it != inps.end(); ++it) {
            outp = intersection(outp, *it);
        }
        return outp;
    }

}