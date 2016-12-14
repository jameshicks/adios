#ifndef SETOPS_HPP
#define SETOPS_HPP

#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <iterator>
#include <initializer_list>


// Shortcut functions for C++'s kind of clumsy set operations in <algorithms> 
namespace setops {
    using std::set;
    using std::unordered_set;
    using std::vector;

    vector<int> union_(const vector<int>& a, const vector<int>& b);
    vector<int> intersection(const vector<int>& a, const vector<int>& b);
    vector<int> symmetric_difference(const vector<int>& a, const vector<int>& b);
    vector<int> difference(const vector<int>& a, const vector<int>& b);
    vector<int> multi_union(std::initializer_list<vector<int>> inps);
    vector<int> multi_intersection(std::initializer_list<vector<int>> inps);

    set<int> union_(const set<int>& a, const set<int>& b);
    set<int> intersection(const set<int>& a, const set<int>& b);
    set<int> symmetric_difference(const set<int>& a, const set<int>& b);
    set<int> difference(const set<int>& a, const set<int>& b);
    set<int> multi_union(std::initializer_list<set<int>> inps);
    set<int> multi_intersection(std::initializer_list<set<int>> inps);
    unordered_set<int> union_(const unordered_set<int>& a, const unordered_set<int>& b);
    unordered_set<int> intersection(const unordered_set<int>& a, const unordered_set<int>& b);
    unordered_set<int> symmetric_difference(const unordered_set<int>& a, const unordered_set<int>& b);
    unordered_set<int> difference(const unordered_set<int>& a, const unordered_set<int>& b);
    unordered_set<int> multi_union(std::initializer_list<unordered_set<int>> inps);
    unordered_set<int> multi_intersection(std::initializer_list<unordered_set<int>> inps);
}

#endif
