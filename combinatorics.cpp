#include "combinatorics.hpp"

namespace combinatorics {

long nCk(long n, long k) {
    double out = 1;
    while (k > 0) {
        out *= n/((double)k); // int division will give wrong answers
        --n;
        --k;
    } 
    return (long)out; // result is always an int
}

long largestV(long a, long b, long x) {
    long v = a - 1;

    while (nCk(v,b) > x) v--;
    return v;
}

std::vector<long> combination_at_index(long i, long n, long k) {
    std::vector<long> pos(k);    

    long a = n;
    long b = k;
    long x = nCk(n, k) - i - 1;

    for (long i = 0; i < k; ++i) {
        pos[i] = largestV(a,b,x);
        x -= nCk(pos[i], b);
        a = pos[i];
        b--;

        pos[i] = n - 1 - pos[i];
    }

    return pos;
}

}