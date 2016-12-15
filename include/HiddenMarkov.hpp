#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>

#include "Linalg.hpp"
using Linalg::Matrix;
using Linalg::Vector;
class GenotypeHMM
{
public:
    int nstates;
    std::vector<int> observations;
    std::vector<Matrix> emission_matrices;
    Matrix transition_matrix;

    GenotypeHMM(const std::vector<int>& obs,
                const std::vector<Matrix>& emission,
                const Matrix& transition);
    inline std::vector<int> decode(bool use_posteriori)
    {
        return use_posteriori ? viterbi() : forwards_backwards();
    }
    std::vector<int> viterbi(void) const;
    std::vector<int> forwards_backwards(void) const;
};

