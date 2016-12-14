#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <vector>
#include <utility>
#include <math.h>


namespace combinatorics
{
using std::pair;

int nCk(int n, int k);

template <typename T>
std::vector<pair<T, T>> pair_combinations(const std::vector<T>& v)
{
    int n = v.size();
    std::vector<pair<T, T>> pairs(nCk(n, 2));

    int index = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            pairs[index] = std::make_pair(v[i], v[j]);
            index++;
        }
    }

    return pairs;
}

}


#endif