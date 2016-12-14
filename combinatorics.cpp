#include "combinatorics.hpp"

namespace combinatorics {

int nCk(int n, int k) {
    double out = 1;
    while (k > 0) {
        out *= n/((double)k); // int division will give wrong answers
        --n;
        --k;
    } 
    return (int)out; // result is always an int
}



}